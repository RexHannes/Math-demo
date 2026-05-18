#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct PrimeInfo {
    int p;
};

struct DenominatorPrimeData {
    int valuation;
    int reduced_inverse;
};

struct CandidateSample {
    vector<int> S;
    double sum_float;
    uint64_t kill_mask;
};

struct ProgressState {
    uint64_t nodes = 0;
    uint64_t near_count = 0;
    vector<int> first_uncovered;
};

static vector<int> primes_upto(int n) {
    if (n < 2) return {};
    vector<bool> sieve(n + 1, true);
    sieve[0] = sieve[1] = false;
    for (int i = 2; i * i <= n; ++i) {
        if (!sieve[i]) continue;
        for (int j = i * i; j <= n; j += i) sieve[j] = false;
    }
    vector<int> primes;
    for (int i = 2; i <= n; ++i) if (sieve[i]) primes.push_back(i);
    return primes;
}

static int mod_pow(int base, int exp, int mod) {
    long long result = 1;
    long long value = base % mod;
    while (exp > 0) {
        if (exp & 1) result = (result * value) % mod;
        value = (value * value) % mod;
        exp >>= 1;
    }
    return static_cast<int>(result);
}

static int mod_inverse_prime(int a, int p) {
    int normalized = a % p;
    if (normalized < 0) normalized += p;
    return mod_pow(normalized, p - 2, p);
}

static vector<double> suffix_harmonic_float(int N) {
    vector<double> suffix(N + 3, 0.0);
    for (int k = N; k >= 2; --k) suffix[k] = suffix[k + 1] + 1.0 / static_cast<double>(k);
    return suffix;
}

static string json_escape(const string& value) {
    string out;
    for (char c : value) {
        if (c == '\\') out += "\\\\";
        else if (c == '"') out += "\\\"";
        else out += c;
    }
    return out;
}

static vector<string> bits_to_labels(uint64_t mask, const vector<string>& labels) {
    vector<string> out;
    for (size_t i = 0; i < labels.size(); ++i) {
        if ((mask >> i) & 1ULL) out.push_back(labels[i]);
    }
    return out;
}

static pair<int, vector<uint64_t>> minimal_hitting_sets(const vector<uint64_t>& unique_masks, int num_tests, int max_solutions = 20) {
    if (unique_masks.empty()) return {0, {0}};

    vector<uint64_t> masks = unique_masks;
    sort(masks.begin(), masks.end(), [](uint64_t left, uint64_t right) {
        int left_bits = __builtin_popcountll(left);
        int right_bits = __builtin_popcountll(right);
        if (left_bits != right_bits) return left_bits < right_bits;
        return left < right;
    });
    masks.erase(unique(masks.begin(), masks.end()), masks.end());

    vector<uint64_t> solutions;
    function<bool(int, int, uint64_t)> search = [&](int start_bit, int remaining, uint64_t cover_mask) {
        if (remaining == 0) {
            for (uint64_t candidate_mask : masks) {
                if ((candidate_mask & cover_mask) == 0) return false;
            }
            solutions.push_back(cover_mask);
            return static_cast<int>(solutions.size()) >= max_solutions;
        }
        for (int bit = start_bit; bit <= num_tests - remaining; ++bit) {
            if (search(bit + 1, remaining - 1, cover_mask | (1ULL << bit))) return true;
        }
        return false;
    };

    for (int size = 1; size <= num_tests; ++size) {
        solutions.clear();
        search(0, size, 0);
        if (!solutions.empty()) return {size, solutions};
    }

    return {num_tests + 1, {}};
}

struct RunConfig {
    int N = 60;
    int P = 60;
    double low = 0.999;
    double high = 1.001;
    string out_path;
    string progress_path;
    uint64_t progress_every = 100000;
    int sample_limit = 5;
};

struct RunResult {
    int N = 0;
    int P = 0;
    double low = 0.0;
    double high = 0.0;
    int test_count = 0;
    vector<string> tests;
    uint64_t nodes = 0;
    uint64_t near_count = 0;
    uint64_t unique_mask_count = 0;
    int minimal_cover_size = 0;
    vector<vector<string>> minimal_covers_sample;
    vector<int> first_uncovered;
    vector<pair<uint64_t, uint64_t>> top_masks;
    vector<CandidateSample> samples;
    double seconds = 0.0;
};

static void write_progress_snapshot(const string& path, const map<string, string>& fields) {
    if (path.empty()) return;
    ofstream out(path);
    out << "{\n";
    bool first = true;
    for (const auto& [key, value] : fields) {
        if (!first) out << ",\n";
        first = false;
        out << "  \"" << json_escape(key) << "\": " << value;
    }
    out << "\n}\n";
}

static RunResult run_cover_hyb(const RunConfig& config) {
    auto started = chrono::steady_clock::now();
    vector<int> primes = primes_upto(config.P);
    vector<string> labels;
    for (int p : primes) labels.push_back(to_string(p));

    vector<double> suffix = suffix_harmonic_float(config.N);
    vector<vector<DenominatorPrimeData>> per_denom(config.N + 1, vector<DenominatorPrimeData>(primes.size()));
    for (int n = 2; n <= config.N; ++n) {
        for (size_t i = 0; i < primes.size(); ++i) {
            int p = primes[i];
            int reduced = n;
            int valuation = 0;
            while (reduced % p == 0) {
                reduced /= p;
                valuation += 1;
            }
            per_denom[n][i] = {valuation, mod_inverse_prime(reduced, p)};
        }
    }

    map<uint64_t, uint64_t> unique_masks;
    vector<CandidateSample> samples;
    ProgressState progress;
    vector<int> current_set;
    vector<int> current_max_v(primes.size(), 0);
    vector<int> current_top_sum(primes.size(), 0);

    auto heartbeat = [&](int last, double current_sum) {
        if (config.progress_every == 0 || progress.nodes % config.progress_every != 0) return;
        double elapsed = chrono::duration<double>(chrono::steady_clock::now() - started).count();
        cout << "progress N=" << config.N
             << ": nodes=" << progress.nodes
             << ", near=" << progress.near_count
             << ", last=" << last
             << ", sum=" << fixed << setprecision(8) << current_sum
             << ", elapsed=" << setprecision(2) << elapsed << "s"
             << endl;

        map<string, string> fields;
        fields["N"] = to_string(config.N);
        fields["P"] = to_string(config.P);
        fields["nodes"] = to_string(progress.nodes);
        fields["near_count"] = to_string(progress.near_count);
        fields["last"] = to_string(last);
        {
            ostringstream out;
            out << fixed << setprecision(8) << current_sum;
            fields["sum_float"] = "\"" + out.str() + "\"";
        }
        {
            ostringstream out;
            out << fixed << setprecision(2) << elapsed;
            fields["elapsed_seconds"] = "\"" + out.str() + "\"";
        }
        if (progress.first_uncovered.empty()) {
            fields["first_uncovered"] = "null";
        } else {
            ostringstream out;
            out << "[";
            for (size_t i = 0; i < progress.first_uncovered.size(); ++i) {
                if (i) out << ", ";
                out << progress.first_uncovered[i];
            }
            out << "]";
            fields["first_uncovered"] = out.str();
        }
        write_progress_snapshot(config.progress_path, fields);
    };

    function<void(int, double, const vector<int>&, const vector<int>&)> visit =
        [&](int last, double current_sum, const vector<int>& max_v, const vector<int>& top_sum) {
            progress.nodes += 1;
            heartbeat(last, current_sum);

            if (current_sum > config.high) return;
            if (last + 1 <= config.N && current_sum + suffix[last + 1] < config.low) return;

            if (config.low <= current_sum && current_sum <= config.high) {
                uint64_t kill_mask = 0;
                for (size_t i = 0; i < primes.size(); ++i) {
                    int p = primes[i];
                    int residue = max_v[i] == 0 ? (top_sum[i] - 1) : top_sum[i];
                    residue %= p;
                    if (residue < 0) residue += p;
                    if (residue != 0) kill_mask |= (1ULL << i);
                }

                unique_masks[kill_mask] += 1;
                progress.near_count += 1;

                if (kill_mask == 0 && progress.first_uncovered.empty()) {
                    progress.first_uncovered = current_set;
                }

                if (static_cast<int>(samples.size()) < config.sample_limit) {
                    samples.push_back({current_set, current_sum, kill_mask});
                }
            }

            for (int next : {last + 1, last + 2}) {
                if (next > config.N) continue;
                vector<int> next_max_v = max_v;
                vector<int> next_top_sum = top_sum;
                for (size_t i = 0; i < primes.size(); ++i) {
                    const auto& datum = per_denom[next][i];
                    if (datum.valuation > next_max_v[i]) {
                        next_max_v[i] = datum.valuation;
                        next_top_sum[i] = datum.reduced_inverse % primes[i];
                    } else if (datum.valuation == next_max_v[i]) {
                        next_top_sum[i] = (next_top_sum[i] + datum.reduced_inverse) % primes[i];
                    }
                }

                current_set.push_back(next);
                visit(next, current_sum + 1.0 / static_cast<double>(next), next_max_v, next_top_sum);
                current_set.pop_back();
            }
        };

    for (int start = 2; start <= config.N; ++start) {
        vector<int> next_max_v = current_max_v;
        vector<int> next_top_sum = current_top_sum;
        for (size_t i = 0; i < primes.size(); ++i) {
            const auto& datum = per_denom[start][i];
            if (datum.valuation > next_max_v[i]) {
                next_max_v[i] = datum.valuation;
                next_top_sum[i] = datum.reduced_inverse % primes[i];
            } else if (datum.valuation == next_max_v[i]) {
                next_top_sum[i] = (next_top_sum[i] + datum.reduced_inverse) % primes[i];
            }
        }
        current_set = {start};
        visit(start, 1.0 / static_cast<double>(start), next_max_v, next_top_sum);
    }

    vector<uint64_t> mask_keys;
    for (const auto& [mask, count] : unique_masks) mask_keys.push_back(mask);
    auto [min_cover_size, cover_masks] = minimal_hitting_sets(mask_keys, static_cast<int>(primes.size()));

    vector<pair<uint64_t, uint64_t>> sorted_masks(unique_masks.begin(), unique_masks.end());
    sort(sorted_masks.begin(), sorted_masks.end(), [](const auto& left, const auto& right) {
        if (left.second != right.second) return left.second > right.second;
        return left.first < right.first;
    });
    if (sorted_masks.size() > 20) sorted_masks.resize(20);

    RunResult result;
    result.N = config.N;
    result.P = config.P;
    result.low = config.low;
    result.high = config.high;
    result.test_count = static_cast<int>(labels.size());
    result.tests = labels;
    result.nodes = progress.nodes;
    result.near_count = progress.near_count;
    result.unique_mask_count = static_cast<uint64_t>(unique_masks.size());
    result.minimal_cover_size = min_cover_size;
    for (uint64_t mask : cover_masks) result.minimal_covers_sample.push_back(bits_to_labels(mask, labels));
    result.first_uncovered = progress.first_uncovered;
    result.top_masks = sorted_masks;
    result.samples = samples;
    result.seconds = chrono::duration<double>(chrono::steady_clock::now() - started).count();
    return result;
}

static void write_result_json(const RunResult& result, const string& out_path) {
    ofstream out(out_path);
    out << "{\n";
    out << "  \"generated_at\": \"" << time(nullptr) << "\",\n";
    out << "  \"runs\": [\n";
    out << "    {\n";
    out << "      \"mode\": \"hyb\",\n";
    out << "      \"N\": " << result.N << ",\n";
    out << "      \"P\": " << result.P << ",\n";
    out << "      \"low\": " << result.low << ",\n";
    out << "      \"high\": " << result.high << ",\n";
    out << "      \"test_count\": " << result.test_count << ",\n";
    out << "      \"tests\": [";
    for (size_t i = 0; i < result.tests.size(); ++i) {
        if (i) out << ", ";
        out << '"' << result.tests[i] << '"';
    }
    out << "],\n";
    out << "      \"nodes\": " << result.nodes << ",\n";
    out << "      \"near_count\": " << result.near_count << ",\n";
    out << "      \"exact_hit_count\": null,\n";
    out << "      \"exact_hits_sample\": [],\n";
    out << "      \"unique_mask_count\": " << result.unique_mask_count << ",\n";
    out << "      \"minimal_cover_size\": " << result.minimal_cover_size << ",\n";
    out << "      \"minimal_covers_sample\": [";
    for (size_t i = 0; i < result.minimal_covers_sample.size(); ++i) {
        if (i) out << ", ";
        out << "[";
        for (size_t j = 0; j < result.minimal_covers_sample[i].size(); ++j) {
            if (j) out << ", ";
            out << '"' << result.minimal_covers_sample[i][j] << '"';
        }
        out << "]";
    }
    out << "],\n";
    out << "      \"first_uncovered\": ";
    if (result.first_uncovered.empty()) {
        out << "null";
    } else {
        out << "[";
        for (size_t i = 0; i < result.first_uncovered.size(); ++i) {
            if (i) out << ", ";
            out << result.first_uncovered[i];
        }
        out << "]";
    }
    out << ",\n";
    out << "      \"top_masks\": [\n";
    for (size_t i = 0; i < result.top_masks.size(); ++i) {
        const auto& [mask, count] = result.top_masks[i];
        out << "        {\"count\": " << count << ", \"kills\": [";
        auto labels = bits_to_labels(mask, result.tests);
        for (size_t j = 0; j < labels.size(); ++j) {
            if (j) out << ", ";
            out << '"' << labels[j] << '"';
        }
        out << "], \"mask\": " << mask << "}";
        if (i + 1 < result.top_masks.size()) out << ",";
        out << "\n";
    }
    out << "      ],\n";
    out << "      \"sample_candidates\": [\n";
    for (size_t i = 0; i < result.samples.size(); ++i) {
        const auto& sample = result.samples[i];
        out << "        {\"S\": [";
        for (size_t j = 0; j < sample.S.size(); ++j) {
            if (j) out << ", ";
            out << sample.S[j];
        }
        out << "], \"sum_float\": " << fixed << setprecision(10) << sample.sum_float
            << ", \"kill_mask\": " << sample.kill_mask << ", \"kills\": [";
        auto labels = bits_to_labels(sample.kill_mask, result.tests);
        for (size_t j = 0; j < labels.size(); ++j) {
            if (j) out << ", ";
            out << '"' << labels[j] << '"';
        }
        out << "]}";
        if (i + 1 < result.samples.size()) out << ",";
        out << "\n";
    }
    out << "      ],\n";
    out << "      \"seconds\": " << result.seconds << "\n";
    out << "    }\n";
    out << "  ]\n";
    out << "}\n";
}

int main(int argc, char** argv) {
    RunConfig config;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        auto next = [&]() { return string(argv[++i]); };
        if (arg == "--N") config.N = stoi(next());
        else if (arg == "--P") config.P = stoi(next());
        else if (arg == "--low") config.low = stod(next());
        else if (arg == "--high") config.high = stod(next());
        else if (arg == "--out") config.out_path = next();
        else if (arg == "--progress-out") config.progress_path = next();
        else if (arg == "--progress-every") config.progress_every = stoull(next());
    }

    if (config.out_path.empty()) {
        cerr << "--out is required\n";
        return 1;
    }

    cout << "Running mode=hyb, N=" << config.N
         << ", P=" << config.P
         << ", window=[" << config.low << "," << config.high << "]"
         << endl;
    RunResult result = run_cover_hyb(config);
    write_result_json(result, config.out_path);

    map<string, string> fields;
    fields["completed"] = "true";
    fields["N"] = to_string(result.N);
    fields["P"] = to_string(result.P);
    fields["nodes"] = to_string(result.nodes);
    fields["near_count"] = to_string(result.near_count);
    fields["unique_mask_count"] = to_string(result.unique_mask_count);
    fields["minimal_cover_size"] = to_string(result.minimal_cover_size);
    fields["seconds"] = "\"" + to_string(result.seconds) + "\"";
    write_progress_snapshot(config.progress_path, fields);

    cout << "{\n"
         << "  \"mode\": \"hyb\",\n"
         << "  \"N\": " << result.N << ",\n"
         << "  \"P\": " << result.P << ",\n"
         << "  \"nodes\": " << result.nodes << ",\n"
         << "  \"near_count\": " << result.near_count << ",\n"
         << "  \"unique_mask_count\": " << result.unique_mask_count << ",\n"
         << "  \"minimal_cover_size\": " << result.minimal_cover_size << ",\n"
         << "  \"first_uncovered\": " << (result.first_uncovered.empty() ? "null" : "\"present\"") << ",\n"
         << "  \"seconds\": " << result.seconds << "\n"
         << "}\n";
    cerr << "Wrote " << config.out_path << "\n";
    return 0;
}
