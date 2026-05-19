#include <algorithm>
#include <array>
#include <chrono>
#include <csignal>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
using i128 = __int128_t;

static volatile sig_atomic_t g_stop_requested = 0;

static void handle_stop_signal(int) {
    g_stop_requested = 1;
}

static bool stop_requested() {
    return g_stop_requested != 0;
}

static void install_signal_handlers() {
    signal(SIGINT, handle_stop_signal);
    signal(SIGTERM, handle_stop_signal);
}

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

struct MonitoredPrimeSnapshot {
    int prime = 0;
    int max_valuation = 0;
    int residue_sum_mod_p = 0;
    vector<int> top_denominators;
};

struct AnomalyRecord {
    uint64_t candidate_id = 0;
    vector<string> escaped_backbones;
    vector<int> S;
    double sum_float = 0.0;
    string exact_sum;
    string sum_minus_one;
    uint64_t kill_mask = 0;
    vector<string> kill_primes;
    vector<MonitoredPrimeSnapshot> snapshots;
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

static string join_strings(const vector<string>& values, const string& separator) {
    ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << separator;
        out << values[i];
    }
    return out.str();
}

static string join_ints(const vector<int>& values, const string& separator) {
    ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << separator;
        out << values[i];
    }
    return out.str();
}

static string csv_escape(const string& value) {
    bool needs_quotes = value.find_first_of(",\"\n") != string::npos;
    if (!needs_quotes) return value;
    string escaped = "\"";
    for (char c : value) {
        if (c == '"') escaped += "\"\"";
        else escaped += c;
    }
    escaped += "\"";
    return escaped;
}

static i128 i128_abs(i128 value) {
    return value < 0 ? -value : value;
}

static i128 i128_gcd(i128 left, i128 right) {
    left = i128_abs(left);
    right = i128_abs(right);
    while (right != 0) {
        i128 next = left % right;
        left = right;
        right = next;
    }
    return left;
}

static i128 i128_pow_int(int base, int exp) {
    i128 value = 1;
    for (int i = 0; i < exp; ++i) value *= base;
    return value;
}

static string i128_to_string(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    i128 current = negative ? -value : value;
    string out;
    while (current > 0) {
        int digit = static_cast<int>(current % 10);
        out.push_back(static_cast<char>('0' + digit));
        current /= 10;
    }
    if (negative) out.push_back('-');
    reverse(out.begin(), out.end());
    return out;
}

static pair<string, string> exact_sum_strings(const vector<int>& denominators) {
    int max_denominator = 0;
    for (int denominator : denominators) max_denominator = max(max_denominator, denominator);
    vector<int> primes = primes_upto(max_denominator);
    vector<int> max_exponents(primes.size(), 0);
    for (int denominator : denominators) {
        int value = denominator;
        for (size_t i = 0; i < primes.size(); ++i) {
            int prime = primes[i];
            int exponent = 0;
            while (value % prime == 0) {
                value /= prime;
                exponent += 1;
            }
            if (exponent > max_exponents[i]) max_exponents[i] = exponent;
            if (value == 1) break;
        }
    }

    i128 lcm = 1;
    for (size_t i = 0; i < primes.size(); ++i) {
        if (max_exponents[i] > 0) lcm *= i128_pow_int(primes[i], max_exponents[i]);
    }

    i128 sum_numerator = 0;
    for (int denominator : denominators) sum_numerator += lcm / denominator;
    i128 delta_numerator = sum_numerator - lcm;

    i128 sum_gcd = i128_gcd(sum_numerator, lcm);
    i128 delta_gcd = i128_gcd(delta_numerator, lcm);
    return {
        i128_to_string(sum_numerator / sum_gcd) + "/" + i128_to_string(lcm / sum_gcd),
        i128_to_string(delta_numerator / delta_gcd) + "/" + i128_to_string(lcm / delta_gcd),
    };
}

static int prime_valuation(int denominator, int prime) {
    int value = denominator;
    int valuation = 0;
    while (value % prime == 0) {
        value /= prime;
        valuation += 1;
    }
    return valuation;
}

static MonitoredPrimeSnapshot build_prime_snapshot(const vector<int>& denominators, int prime) {
    MonitoredPrimeSnapshot snapshot;
    snapshot.prime = prime;
    for (int denominator : denominators) {
        snapshot.max_valuation = max(snapshot.max_valuation, prime_valuation(denominator, prime));
    }
    int residue_sum = 0;
    for (int denominator : denominators) {
        int valuation = prime_valuation(denominator, prime);
        if (valuation != snapshot.max_valuation) continue;
        snapshot.top_denominators.push_back(denominator);
        int reduced = denominator;
        for (int i = 0; i < valuation; ++i) reduced /= prime;
        residue_sum = (residue_sum + mod_inverse_prime(reduced, prime)) % prime;
    }
    snapshot.residue_sum_mod_p = residue_sum;
    return snapshot;
}

static vector<int> gap_pattern(const vector<int>& denominators) {
    vector<int> gaps;
    for (size_t i = 1; i < denominators.size(); ++i) gaps.push_back(denominators[i] - denominators[i - 1]);
    return gaps;
}

static bool contains_pair(const vector<int>& denominators, int left, int right) {
    bool seen_left = false;
    bool seen_right = false;
    for (int denominator : denominators) {
        if (denominator == left) seen_left = true;
        if (denominator == right) seen_right = true;
    }
    return seen_left && seen_right;
}

static bool contains_adjacent_factor_pair(const vector<int>& denominators) {
    for (size_t i = 1; i < denominators.size(); ++i) {
        int left = denominators[i - 1];
        int right = denominators[i];
        if (right % left == 0) return true;
    }
    return false;
}

struct BackboneCheck {
    int left = 0;
    int right = 0;
    int left_index = -1;
    int right_index = -1;
};

static string backbone_label(const BackboneCheck& backbone) {
    return to_string(backbone.left) + "+" + to_string(backbone.right);
}

static vector<BackboneCheck> default_backbones(const vector<int>& primes) {
    vector<BackboneCheck> out;
    for (const auto& pair : array<pair<int, int>, 2>{pair<int, int>{2, 31}, pair<int, int>{19, 37}}) {
        BackboneCheck backbone;
        backbone.left = pair.first;
        backbone.right = pair.second;
        for (size_t i = 0; i < primes.size(); ++i) {
            if (primes[i] == backbone.left) backbone.left_index = static_cast<int>(i);
            if (primes[i] == backbone.right) backbone.right_index = static_cast<int>(i);
        }
        out.push_back(backbone);
    }
    return out;
}

static vector<string> escaping_backbones(uint64_t kill_mask, const vector<BackboneCheck>& backbones) {
    vector<string> out;
    for (const auto& backbone : backbones) {
        if (backbone.left_index < 0 || backbone.right_index < 0) continue;
        bool left_hit = ((kill_mask >> backbone.left_index) & 1ULL) != 0;
        bool right_hit = ((kill_mask >> backbone.right_index) & 1ULL) != 0;
        if (!left_hit && !right_hit) out.push_back(backbone_label(backbone));
    }
    return out;
}

static string covering_third_primes_string(const vector<string>& escaped_backbones, const vector<string>& kill_primes) {
    vector<string> mappings;
    string joined_primes = kill_primes.empty() ? "none" : join_strings(kill_primes, ";");
    for (const auto& backbone : escaped_backbones) mappings.push_back(backbone + "->" + joined_primes);
    return join_strings(mappings, "|");
}

static const vector<int>& monitored_primes() {
    static const vector<int> primes{2, 13, 17, 19, 31, 37, 41, 47, 53, 59, 61};
    return primes;
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
    int start_min = 2;
    int start_max = 0;
    string out_path;
    string progress_path;
    string anomalies_csv_path;
    string anomaly_summary_path;
    string chunk_status_path;
    string chunk_id;
    int num_chunks = 0;
    uint64_t progress_every = 10000000;
    int sample_limit = 5;
};

struct RunResult {
    bool completed = true;
    bool interrupted = false;
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
    vector<pair<uint64_t, uint64_t>> unique_masks;
    vector<CandidateSample> samples;
    vector<AnomalyRecord> anomalies;
    double seconds = 0.0;
};

static void write_json_atomically(const string& path, const function<void(ostream&)>& writer) {
    if (path.empty()) return;
    string tmp_path = path + ".tmp";
    ofstream out(tmp_path);
    writer(out);
    out.close();
    rename(tmp_path.c_str(), path.c_str());
}

static void write_progress_snapshot(const string& path, const map<string, string>& fields) {
    if (path.empty()) return;
    write_json_atomically(path, [&](ostream& out) {
        out << "{\n";
        bool first = true;
        for (const auto& [key, value] : fields) {
            if (!first) out << ",\n";
            first = false;
            out << "  \"" << json_escape(key) << "\": " << value;
        }
        out << "\n}\n";
    });
}

static void write_partial_result_json(const RunConfig& config,
                                      const vector<string>& tests,
                                      const ProgressState& progress,
                                      const vector<CandidateSample>& samples,
                                      const unordered_map<uint64_t, uint64_t>& unique_masks,
                                      const string& out_path,
                                      double seconds,
                                      bool interrupted,
                                      int checkpoint_last,
                                      double checkpoint_sum) {
    vector<pair<uint64_t, uint64_t>> sorted_unique_masks(unique_masks.begin(), unique_masks.end());
    sort(sorted_unique_masks.begin(), sorted_unique_masks.end());

    write_json_atomically(out_path, [&](ostream& out) {
        out << "{\n";
        out << "  \"generated_at\": \"" << time(nullptr) << "\",\n";
        out << "  \"completed\": false,\n";
        out << "  \"interrupted\": " << (interrupted ? "true" : "false") << ",\n";
        out << "  \"runs\": [\n";
        out << "    {\n";
        out << "      \"mode\": \"hyb\",\n";
        out << "      \"N\": " << config.N << ",\n";
        out << "      \"P\": " << config.P << ",\n";
        out << "      \"low\": " << config.low << ",\n";
        out << "      \"high\": " << config.high << ",\n";
        out << "      \"start_min\": " << config.start_min << ",\n";
        out << "      \"start_max\": " << (config.start_max > 0 ? config.start_max : config.N) << ",\n";
        out << "      \"test_count\": " << tests.size() << ",\n";
        out << "      \"tests\": [";
        for (size_t i = 0; i < tests.size(); ++i) {
            if (i) out << ", ";
            out << '"' << tests[i] << '"';
        }
        out << "],\n";
        out << "      \"nodes\": " << progress.nodes << ",\n";
        out << "      \"near_count\": " << progress.near_count << ",\n";
        out << "      \"exact_hit_count\": null,\n";
        out << "      \"exact_hits_sample\": [],\n";
        out << "      \"unique_mask_count\": " << sorted_unique_masks.size() << ",\n";
        out << "      \"unique_masks\": {";
        for (size_t i = 0; i < sorted_unique_masks.size(); ++i) {
            const auto& [mask, count] = sorted_unique_masks[i];
            if (i) out << ", ";
            out << '"' << mask << "\": " << count;
        }
        out << "},\n";
        out << "      \"minimal_cover_size\": null,\n";
        out << "      \"minimal_covers_sample\": [],\n";
        out << "      \"first_uncovered\": ";
        if (progress.first_uncovered.empty()) {
            out << "null";
        } else {
            out << "[";
            for (size_t i = 0; i < progress.first_uncovered.size(); ++i) {
                if (i) out << ", ";
                out << progress.first_uncovered[i];
            }
            out << "]";
        }
        out << ",\n";
        out << "      \"top_masks\": [],\n";
        out << "      \"sample_candidates\": [\n";
        for (size_t i = 0; i < samples.size(); ++i) {
            const auto& sample = samples[i];
            out << "        {\"S\": [";
            for (size_t j = 0; j < sample.S.size(); ++j) {
                if (j) out << ", ";
                out << sample.S[j];
            }
            out << "], \"sum_float\": " << fixed << setprecision(10) << sample.sum_float
                << ", \"kill_mask\": " << sample.kill_mask << ", \"kills\": [";
            auto labels = bits_to_labels(sample.kill_mask, tests);
            for (size_t j = 0; j < labels.size(); ++j) {
                if (j) out << ", ";
                out << '"' << labels[j] << '"';
            }
            out << "]}";
            if (i + 1 < samples.size()) out << ",";
            out << "\n";
        }
        out << "      ],\n";
        out << "      \"checkpoint_last\": " << checkpoint_last << ",\n";
        out << "      \"checkpoint_sum_float\": " << fixed << setprecision(10) << checkpoint_sum << ",\n";
        out << "      \"seconds\": " << seconds << "\n";
        out << "    }\n";
        out << "  ]\n";
        out << "}\n";
    });
}

static AnomalyRecord build_anomaly_record(uint64_t candidate_id,
                                          const vector<int>& denominators,
                                          double sum_float,
                                          uint64_t kill_mask,
                                          const vector<string>& kill_primes,
                                          const vector<string>& escaped_backbones) {
    AnomalyRecord record;
    record.candidate_id = candidate_id;
    record.escaped_backbones = escaped_backbones;
    record.S = denominators;
    record.sum_float = sum_float;
    record.kill_mask = kill_mask;
    record.kill_primes = kill_primes;
    auto fractions = exact_sum_strings(denominators);
    record.exact_sum = fractions.first;
    record.sum_minus_one = fractions.second;

    for (int prime : monitored_primes()) {
        if (prime <= 1) continue;
        record.snapshots.push_back(build_prime_snapshot(denominators, prime));
    }
    return record;
}

static void write_anomaly_csv_header(ostream& out) {
    out << "candidate_id,escaped_backbone,denominators,min_denominator,max_denominator,length,gap_pattern,exact_sum,sum_minus_1,kill_mask,kill_primes,kill_count,first_killing_prime,covering_third_primes,contains_37_38,contains_61_62,contains_adjacent_factor_pair";
    for (int prime : monitored_primes()) {
        out << ",valuation_profile_p" << prime << ",top_layer_residue_sum_p" << prime;
    }
    out << "\n";
}

static void write_anomaly_csv_record(ostream& out, const AnomalyRecord& record) {
    vector<int> gaps = gap_pattern(record.S);
    vector<string> escaped = record.escaped_backbones;
    vector<string> kill_primes = record.kill_primes;
    string denominators = join_ints(record.S, ";");
    string gap_text = join_ints(gaps, ";");
    string escaped_text = join_strings(escaped, "|");
    string kill_text = join_strings(kill_primes, ";");
    string third_primes = covering_third_primes_string(escaped, kill_primes);
    string first_killing_prime = kill_primes.empty() ? "" : kill_primes.front();

    out << record.candidate_id
        << "," << csv_escape(escaped_text)
        << "," << csv_escape(denominators)
        << "," << record.S.front()
        << "," << record.S.back()
        << "," << record.S.size()
        << "," << csv_escape(gap_text)
        << "," << csv_escape(record.exact_sum)
        << "," << csv_escape(record.sum_minus_one)
        << "," << record.kill_mask
        << "," << csv_escape(kill_text)
        << "," << kill_primes.size()
        << "," << csv_escape(first_killing_prime)
        << "," << csv_escape(third_primes)
        << "," << (contains_pair(record.S, 37, 38) ? "true" : "false")
        << "," << (contains_pair(record.S, 61, 62) ? "true" : "false")
        << "," << (contains_adjacent_factor_pair(record.S) ? "true" : "false");

    for (const auto& snapshot : record.snapshots) {
        vector<string> profile_parts;
        profile_parts.push_back("max=" + to_string(snapshot.max_valuation));
        profile_parts.push_back("tops=" + join_ints(snapshot.top_denominators, ";"));
        out << "," << csv_escape(join_strings(profile_parts, ";"))
            << "," << snapshot.residue_sum_mod_p;
    }
    out << "\n";
}

static RunResult run_cover_hyb(const RunConfig& config) {
    auto started = chrono::steady_clock::now();
    vector<int> primes = primes_upto(config.P);
    vector<BackboneCheck> anomaly_backbones = default_backbones(primes);
    vector<string> labels;
    for (int p : primes) labels.push_back(to_string(p));

    vector<double> reciprocal(config.N + 1, 0.0);
    for (int n = 2; n <= config.N; ++n) reciprocal[n] = 1.0 / static_cast<double>(n);

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

    unordered_map<uint64_t, uint64_t> unique_masks;
    if (primes.size() < 20) {
        unique_masks.reserve(1ULL << primes.size());
    } else {
        unique_masks.reserve(1ULL << 20);
    }
    vector<CandidateSample> samples;
    vector<AnomalyRecord> anomalies;
    ProgressState progress;
    vector<int> current_set;
    vector<int> current_max_v(primes.size(), 0);
    vector<int> current_top_sum(primes.size(), 0);
    int last_checkpoint_last = 0;
    double last_checkpoint_sum = 0.0;
    ofstream anomaly_csv;
    bool stream_anomalies = !config.anomalies_csv_path.empty();
    if (stream_anomalies) {
        anomaly_csv.open(config.anomalies_csv_path);
        if (!anomaly_csv) {
            cerr << "Could not open anomaly CSV for streaming: " << config.anomalies_csv_path << "\n";
            stream_anomalies = false;
        } else {
            write_anomaly_csv_header(anomaly_csv);
            anomaly_csv.flush();
        }
    }

    auto checkpoint_partial = [&](int last, double current_sum, bool interrupted) {
        last_checkpoint_last = last;
        last_checkpoint_sum = current_sum;
        double elapsed = chrono::duration<double>(chrono::steady_clock::now() - started).count();
        write_partial_result_json(
            config,
            labels,
            progress,
            samples,
            unique_masks,
            config.out_path,
            elapsed,
            interrupted,
            last,
            current_sum
        );
    };

    checkpoint_partial(0, 0.0, false);

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
        checkpoint_partial(last, current_sum, stop_requested());
    };

    function<void(int, double, const vector<int>&, const vector<int>&)> visit =
        [&](int last, double current_sum, const vector<int>& max_v, const vector<int>& top_sum) {
            if (stop_requested()) return;
            progress.nodes += 1;
            heartbeat(last, current_sum);

            if (stop_requested()) return;
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
                uint64_t candidate_id = progress.near_count + 1;
                progress.near_count += 1;

                if (kill_mask == 0 && progress.first_uncovered.empty()) {
                    progress.first_uncovered = current_set;
                }

                if (static_cast<int>(samples.size()) < config.sample_limit) {
                    samples.push_back({current_set, current_sum, kill_mask});
                }

                if (!config.anomalies_csv_path.empty()) {
                    vector<string> escaped = escaping_backbones(kill_mask, anomaly_backbones);
                    if (!escaped.empty()) {
                        vector<string> kill_primes = bits_to_labels(kill_mask, labels);
                        AnomalyRecord record = build_anomaly_record(
                            candidate_id,
                            current_set,
                            current_sum,
                            kill_mask,
                            kill_primes,
                            escaped
                        );
                        if (stream_anomalies) {
                            write_anomaly_csv_record(anomaly_csv, record);
                            anomaly_csv.flush();
                        }
                        anomalies.push_back(record);
                    }
                }
            }

            for (int next : {last + 1, last + 2}) {
                if (stop_requested()) return;
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
                visit(next, current_sum + reciprocal[next], next_max_v, next_top_sum);
                current_set.pop_back();
            }
        };

    int start_min = max(2, config.start_min);
    int start_max = config.start_max > 0 ? min(config.N, config.start_max) : config.N;
    for (int start = start_min; start <= start_max; ++start) {
        if (stop_requested()) break;
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
        visit(start, reciprocal[start], next_max_v, next_top_sum);
    }

    if (stop_requested()) {
        checkpoint_partial(last_checkpoint_last, last_checkpoint_sum, true);
        RunResult result;
        result.completed = false;
        result.interrupted = true;
        result.N = config.N;
        result.P = config.P;
        result.low = config.low;
        result.high = config.high;
        result.test_count = static_cast<int>(labels.size());
        result.tests = labels;
        result.nodes = progress.nodes;
        result.near_count = progress.near_count;
        result.unique_mask_count = static_cast<uint64_t>(unique_masks.size());
        result.unique_masks.assign(unique_masks.begin(), unique_masks.end());
        sort(result.unique_masks.begin(), result.unique_masks.end());
        result.first_uncovered = progress.first_uncovered;
        result.samples = samples;
        result.anomalies = anomalies;
        result.seconds = chrono::duration<double>(chrono::steady_clock::now() - started).count();
        return result;
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
    result.completed = true;
    result.interrupted = false;
    result.N = config.N;
    result.P = config.P;
    result.low = config.low;
    result.high = config.high;
    result.test_count = static_cast<int>(labels.size());
    result.tests = labels;
    result.nodes = progress.nodes;
    result.near_count = progress.near_count;
    result.unique_mask_count = static_cast<uint64_t>(unique_masks.size());
    result.unique_masks.assign(unique_masks.begin(), unique_masks.end());
    sort(result.unique_masks.begin(), result.unique_masks.end());
    result.minimal_cover_size = min_cover_size;
    for (uint64_t mask : cover_masks) result.minimal_covers_sample.push_back(bits_to_labels(mask, labels));
    result.first_uncovered = progress.first_uncovered;
    result.top_masks = sorted_masks;
    result.samples = samples;
    result.anomalies = anomalies;
    result.seconds = chrono::duration<double>(chrono::steady_clock::now() - started).count();
    return result;
}

static void write_result_json(const RunResult& result, const string& out_path) {
    write_json_atomically(out_path, [&](ostream& out) {
        out << "{\n";
        out << "  \"generated_at\": \"" << time(nullptr) << "\",\n";
        out << "  \"completed\": " << (result.completed ? "true" : "false") << ",\n";
        out << "  \"interrupted\": " << (result.interrupted ? "true" : "false") << ",\n";
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
        out << "      \"unique_masks\": {";
        for (size_t i = 0; i < result.unique_masks.size(); ++i) {
            const auto& [mask, count] = result.unique_masks[i];
            if (i) out << ", ";
            out << '"' << mask << "\": " << count;
        }
        out << "},\n";
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
    });
}

static void write_anomaly_outputs(const RunResult& result,
                                  const string& csv_path,
                                  const string& summary_path) {
    if (csv_path.empty()) return;

    if (summary_path.empty()) return;
    write_json_atomically(summary_path, [&](ostream& out) {
        map<string, int> backbone_counts;
        int contains_37_38 = 0;
        int contains_61_62 = 0;
        int adjacent_factor_pair = 0;

        for (const auto& record : result.anomalies) {
            for (const auto& backbone : record.escaped_backbones) backbone_counts[backbone] += 1;
            if (contains_pair(record.S, 37, 38)) contains_37_38 += 1;
            if (contains_pair(record.S, 61, 62)) contains_61_62 += 1;
            if (contains_adjacent_factor_pair(record.S)) adjacent_factor_pair += 1;
        }

        out << "# Anomaly Candidate Summary\n\n";
        out << "- completed: `" << (result.completed ? "true" : "false") << "`\n";
        out << "- interrupted: `" << (result.interrupted ? "true" : "false") << "`\n";
        out << "- anomaly_rows: `" << result.anomalies.size() << "`\n";
        out << "- near_count: `" << result.near_count << "`\n";
        out << "- output_csv: `" << csv_path << "`\n\n";
        out << "## Escaped Backbones\n\n";
        for (const auto& [backbone, count] : backbone_counts) {
            out << "- `" << backbone << "`: `" << count << "`\n";
        }
        if (backbone_counts.empty()) out << "- none\n";
        out << "\n## Structural Flags\n\n";
        out << "- contains_37_38: `" << contains_37_38 << "`\n";
        out << "- contains_61_62: `" << contains_61_62 << "`\n";
        out << "- contains_adjacent_factor_pair: `" << adjacent_factor_pair << "`\n";
        out << "\n## Sample Rows\n\n";
        size_t sample_count = min<size_t>(5, result.anomalies.size());
        for (size_t i = 0; i < sample_count; ++i) {
            const auto& record = result.anomalies[i];
            out << "- candidate `" << record.candidate_id << "` "
                << "backbones=`" << join_strings(record.escaped_backbones, "|") << "` "
                << "denominators=`" << join_ints(record.S, ",") << "` "
                << "kill_primes=`" << join_strings(record.kill_primes, ",") << "` "
                << "sum_minus_1=`" << record.sum_minus_one << "`\n";
        }
    });
}

static void write_chunk_status(const RunConfig& config, const RunResult& result) {
    if (config.chunk_status_path.empty()) return;
    map<string, int> backbone_counts;
    for (const auto& record : result.anomalies) {
        for (const auto& backbone : record.escaped_backbones) backbone_counts[backbone] += 1;
    }

    write_json_atomically(config.chunk_status_path, [&](ostream& out) {
        out << "{\n";
        out << "  \"N\": " << result.N << ",\n";
        out << "  \"P\": " << result.P << ",\n";
        out << "  \"low\": " << result.low << ",\n";
        out << "  \"high\": " << result.high << ",\n";
        out << "  \"chunk_id\": \"" << json_escape(config.chunk_id.empty() ? "S" + to_string(config.start_min) + "_" + to_string(config.start_max > 0 ? config.start_max : config.N) : config.chunk_id) << "\",\n";
        out << "  \"num_chunks\": " << config.num_chunks << ",\n";
        out << "  \"start_min\": " << config.start_min << ",\n";
        out << "  \"start_max\": " << (config.start_max > 0 ? config.start_max : config.N) << ",\n";
        out << "  \"status\": \"" << (result.completed ? "complete" : "incomplete") << "\",\n";
        out << "  \"completed\": " << (result.completed ? "true" : "false") << ",\n";
        out << "  \"interrupted\": " << (result.interrupted ? "true" : "false") << ",\n";
        out << "  \"nodes\": " << result.nodes << ",\n";
        out << "  \"candidates_checked\": " << result.near_count << ",\n";
        out << "  \"exceptions_found\": " << result.anomalies.size() << ",\n";
        out << "  \"csv_path\": \"" << json_escape(config.anomalies_csv_path) << "\",\n";
        out << "  \"pairs\": [\"2+31\", \"19+37\"],\n";
        out << "  \"backbone_counts\": {";
        bool first = true;
        for (const auto& [backbone, count] : backbone_counts) {
            if (!first) out << ", ";
            first = false;
            out << "\"" << json_escape(backbone) << "\": " << count;
        }
        out << "}\n";
        out << "}\n";
    });
}

int main(int argc, char** argv) {
    install_signal_handlers();
    RunConfig config;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        auto next = [&]() { return string(argv[++i]); };
        if (arg == "--N") config.N = stoi(next());
        else if (arg == "--P") config.P = stoi(next());
        else if (arg == "--low") config.low = stod(next());
        else if (arg == "--high") config.high = stod(next());
        else if (arg == "--start-min") config.start_min = stoi(next());
        else if (arg == "--start-max") config.start_max = stoi(next());
        else if (arg == "--out") config.out_path = next();
        else if (arg == "--progress-out") config.progress_path = next();
        else if (arg == "--progress-every") config.progress_every = stoull(next());
        else if (arg == "--dump-anomalies") config.anomalies_csv_path = next();
        else if (arg == "--anomaly-summary") config.anomaly_summary_path = next();
        else if (arg == "--chunk-status") config.chunk_status_path = next();
        else if (arg == "--chunk-id") config.chunk_id = next();
        else if (arg == "--num-chunks") config.num_chunks = stoi(next());
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
    if (result.completed) {
        write_result_json(result, config.out_path);
    }
    if (!config.anomalies_csv_path.empty()) {
        write_anomaly_outputs(result, config.anomalies_csv_path, config.anomaly_summary_path);
    }
    write_chunk_status(config, result);

    map<string, string> fields;
    fields["completed"] = result.completed ? "true" : "false";
    fields["interrupted"] = result.interrupted ? "true" : "false";
    fields["N"] = to_string(result.N);
    fields["P"] = to_string(result.P);
    fields["nodes"] = to_string(result.nodes);
    fields["near_count"] = to_string(result.near_count);
    fields["unique_mask_count"] = to_string(result.unique_mask_count);
    fields["minimal_cover_size"] = result.completed ? to_string(result.minimal_cover_size) : "null";
    fields["seconds"] = "\"" + to_string(result.seconds) + "\"";
    write_progress_snapshot(config.progress_path, fields);

    cout << "{\n"
         << "  \"completed\": " << (result.completed ? "true" : "false") << ",\n"
         << "  \"interrupted\": " << (result.interrupted ? "true" : "false") << ",\n"
         << "  \"mode\": \"hyb\",\n"
         << "  \"N\": " << result.N << ",\n"
         << "  \"P\": " << result.P << ",\n"
         << "  \"nodes\": " << result.nodes << ",\n"
         << "  \"near_count\": " << result.near_count << ",\n"
         << "  \"unique_mask_count\": " << result.unique_mask_count << ",\n"
         << "  \"anomaly_count\": " << result.anomalies.size() << ",\n"
         << "  \"minimal_cover_size\": " << (result.completed ? to_string(result.minimal_cover_size) : "null") << ",\n"
         << "  \"first_uncovered\": " << (result.first_uncovered.empty() ? "null" : "\"present\"") << ",\n"
         << "  \"seconds\": " << result.seconds << "\n"
         << "}\n";
    cerr << "Wrote " << config.out_path << "\n";
    if (!config.anomalies_csv_path.empty()) cerr << "Wrote " << config.anomalies_csv_path << "\n";
    return 0;
}
