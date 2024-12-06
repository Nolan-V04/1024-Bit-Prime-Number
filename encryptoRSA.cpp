#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <chrono>

using namespace std;

struct BigInt {
    bool sign;  // true for positive, false for negative
    vector<uint32_t> BigIntDigits;
};

const uint32_t BASE = 1000000000; // Base for storing digits (9 decimal digits)

// Convert hex character to integer
uint32_t hexToInt(char hexChar) {
    if (hexChar >= '0' && hexChar <= '9') {
        return hexChar - '0';
    } else if (hexChar >= 'A' && hexChar <= 'F') {
        return hexChar - 'A' + 10;
    } else if (hexChar >= 'a' && hexChar <= 'f') {
        return hexChar - 'a' + 10;
    } else {
        throw invalid_argument("Invalid hex character");
    }
}

// Multiply BigInt by 16
void multiplyBy16(BigInt& bigInt) {
    uint64_t carry = 0;
    for (auto& digit : bigInt.BigIntDigits) {
        uint64_t product = static_cast<uint64_t>(digit) * 16 + carry;
        digit = product % BASE;
        carry = product / BASE;
    }
    if (carry > 0) {
        bigInt.BigIntDigits.push_back(static_cast<uint32_t>(carry));
    }
}

// Add a single digit to BigInt
void addDigit(BigInt& bigInt, uint32_t digit) {
    uint64_t carry = digit;
    for (auto& element : bigInt.BigIntDigits) {
        uint64_t sum = static_cast<uint64_t>(element) + carry;
        element = sum % BASE;
        carry = sum / BASE;
        if (carry == 0) return;
    }
    if (carry > 0) {
        bigInt.BigIntDigits.push_back(static_cast<uint32_t>(carry));
    }
}

// Convert hex string to BigInt
BigInt hexToBigInt(const string& hex) {
    BigInt result = {true, {0}};
    for (char hexChar : hex) {
        multiplyBy16(result);
        addDigit(result, hexToInt(hexChar));
    }
    return result;
}

// Print BigInt
void printBigInt(const BigInt& bigInt) {
    if (bigInt.BigIntDigits.empty()) {
        cout << "0\n";
        return;
    }
    if (!bigInt.sign) cout << "-";

    // Print the most significant digit without leading zeros
    cout << bigInt.BigIntDigits.back();

    // Print the remaining digits with leading zeros
    for (int i = bigInt.BigIntDigits.size() - 2; i >= 0; --i) {
        cout << setw(9) << setfill('0') << bigInt.BigIntDigits[i];
    }
    cout << endl;
}

// Compare two BigInts: -1 if a < b, 0 if a == b, 1 if a > b
int compareBigInts(const BigInt& a, const BigInt& b) {
    // Xét dấu trước
    if (a.sign != b.sign) return a.sign ? 1 : -1;

    // Nếu cùng dấu, so sánh phần giá trị
    int magnitudeComparison = 0;
    if (a.BigIntDigits.size() > b.BigIntDigits.size()) magnitudeComparison = 1;
    else if (a.BigIntDigits.size() < b.BigIntDigits.size()) magnitudeComparison = -1;
    else {
        for (size_t i = a.BigIntDigits.size(); i-- > 0;) {
            if (a.BigIntDigits[i] > b.BigIntDigits[i]) {
                magnitudeComparison = 1;
                break;
            }
            if (a.BigIntDigits[i] < b.BigIntDigits[i]) {
                magnitudeComparison = -1;
                break;
            }
        }
    }

    // Nếu là số âm, đảo kết quả so sánh
    return a.sign ? magnitudeComparison : -magnitudeComparison;
}

BigInt subtractBigInts(const BigInt& a, const BigInt& b);

BigInt addBigInts(const BigInt& a, const BigInt& b) {
    if (a.sign == b.sign) {
        // Nếu cùng dấu, thực hiện cộng
        BigInt result = {a.sign, {}};
        uint64_t carry = 0;
        size_t maxSize = max(a.BigIntDigits.size(), b.BigIntDigits.size());
        result.BigIntDigits.resize(maxSize);

        for (size_t i = 0; i < maxSize; ++i) {
            uint64_t sum = carry;
            if (i < a.BigIntDigits.size()) sum += a.BigIntDigits[i];
            if (i < b.BigIntDigits.size()) sum += b.BigIntDigits[i];
            result.BigIntDigits[i] = sum % BASE;
            carry = sum / BASE;
        }
        if (carry > 0) {
            result.BigIntDigits.push_back(static_cast<uint32_t>(carry));
        }
        return result;
    } else {
        // Nếu khác dấu, chuyển về phép trừ
        if (!a.sign) {
            // a âm, b dương => b - |a|
            return subtractBigInts(b, {true, a.BigIntDigits});
        } else {
            // a dương, b âm => a - |b|
            return subtractBigInts(a, {true, b.BigIntDigits});
        }
    }
}


BigInt subtractBigInts(const BigInt& a, const BigInt& b) {
    if (!a.sign && !b.sign) {
        // Nếu cả hai là số âm, đảo phép trừ: (-a) - (-b) => b - a
        return subtractBigInts(b, a);
    }

    if (a.sign != b.sign) {
        // Nếu khác dấu: a - (-b) => a + b
        return addBigInts(a, {a.sign, b.BigIntDigits});
    }

    // Nếu cùng dấu, kiểm tra giá trị lớn hơn
    if (compareBigInts(a, b) < 0) {
        BigInt result = subtractBigInts(b, a);
        result.sign = !a.sign; // Kết quả đảo dấu
        return result;
    }

    // Thực hiện trừ thông thường
    BigInt result = a;
    int64_t borrow = 0;

    for (size_t i = 0; i < b.BigIntDigits.size() || borrow > 0; ++i) {
        int64_t sub = borrow + (i < b.BigIntDigits.size() ? b.BigIntDigits[i] : 0);
        if (result.BigIntDigits[i] < sub) {
            result.BigIntDigits[i] += BASE - sub;
            borrow = 1;
        } else {
            result.BigIntDigits[i] -= sub;
            borrow = 0;
        }
    }

    // Xóa các chữ số 0 dư thừa
    while (result.BigIntDigits.size() > 1 && result.BigIntDigits.back() == 0) {
        result.BigIntDigits.pop_back();
    }

    return result;
}

BigInt multiplyBigInts(const BigInt& a, const BigInt& b) {
    if (a.BigIntDigits.empty() || b.BigIntDigits.empty()) return {true, {0}};
    BigInt result = {true, vector<uint32_t>(a.BigIntDigits.size() + b.BigIntDigits.size(), 0)};

    // Handle sign
    result.sign = (a.sign == b.sign);

    for (size_t i = 0; i < a.BigIntDigits.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.BigIntDigits.size(); ++j) {
            uint64_t product = static_cast<uint64_t>(a.BigIntDigits[i]) * b.BigIntDigits[j] + result.BigIntDigits[i + j] + carry;
            result.BigIntDigits[i + j] = product % BASE;
            carry = product / BASE;
        }
        if (carry > 0) {
            result.BigIntDigits[i + b.BigIntDigits.size()] += static_cast<uint32_t>(carry);
        }
    }

    // Remove leading zeros
    while (result.BigIntDigits.size() > 1 && result.BigIntDigits.back() == 0) {
        result.BigIntDigits.pop_back();
    }

    return result;
}

pair<BigInt, BigInt> divideBigInts(const BigInt& dividend, const BigInt& divisor) {
    if (divisor.BigIntDigits.empty() || (divisor.BigIntDigits.size() == 1 && divisor.BigIntDigits[0] == 0)) {
        throw invalid_argument("Division by zero");
    }

    if (compareBigInts(dividend, divisor) < 0) {
        return {{true, {0}}, dividend}; // Thương = 0, Dư = dividend
    }

    BigInt quotient = {dividend.sign == divisor.sign, vector<uint32_t>(dividend.BigIntDigits.size(), 0)};
    BigInt remainder = {dividend.sign, dividend.BigIntDigits};

    size_t shift = dividend.BigIntDigits.size() - divisor.BigIntDigits.size();
    BigInt shiftedDivisor = divisor;
    shiftedDivisor.BigIntDigits.insert(shiftedDivisor.BigIntDigits.begin(), shift, 0);

    for (size_t i = shift + 1; i-- > 0;) {
        uint32_t q = 0;
        uint32_t low = 0, high = BASE - 1;

        while (low <= high) {
            uint32_t mid = low + (high - low) / 2;
            BigInt testProduct = multiplyBigInts(shiftedDivisor, {true, {mid}});

            if (compareBigInts(testProduct, remainder) <= 0) {
                q = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }

        quotient.BigIntDigits[i] = q;
        remainder = subtractBigInts(remainder, multiplyBigInts(shiftedDivisor, {true, {q}}));

        if (i > 0) {
            shiftedDivisor.BigIntDigits.erase(shiftedDivisor.BigIntDigits.begin());
        }
    }

    // Xóa các chữ số 0 dư thừa
    while (quotient.BigIntDigits.size() > 1 && quotient.BigIntDigits.back() == 0) {
        quotient.BigIntDigits.pop_back();
    }
    while (remainder.BigIntDigits.size() > 1 && remainder.BigIntDigits.back() == 0) {
        remainder.BigIntDigits.pop_back();
    }

    return {quotient, remainder};
}

// Hàm modBigInts
BigInt modBigInts(const BigInt& a, const BigInt& b) {
    if (b.BigIntDigits.empty() || (b.BigIntDigits.size() == 1 && b.BigIntDigits[0] == 0)) {
        throw invalid_argument("Modulo by zero");
    }

    BigInt remainder = divideBigInts(a, b).second;

    // Nếu phần dư âm, cộng thêm `b`
    if (!remainder.sign) {
        remainder = addBigInts(remainder, b);
    }

    return remainder;
}

BigInt modExp(const BigInt& base, const BigInt& exponent, const BigInt& modulus) {
    if (modulus.BigIntDigits.empty() || (modulus.BigIntDigits.size() == 1 && modulus.BigIntDigits[0] == 0)) {
        throw invalid_argument("Modulo by zero");
    }

    BigInt result = {true, {1}}; // Initialize result as 1
    BigInt baseMod = modBigInts(base, modulus); // Ensure base is within modulus range
    BigInt exp = exponent; // Copy of the exponent for modification

    while (!exp.BigIntDigits.empty() && !(exp.BigIntDigits.size() == 1 && exp.BigIntDigits[0] == 0)) {
        // Check if the current exponent is odd
        if (exp.BigIntDigits[0] % 2 == 1) {
            result = modBigInts(multiplyBigInts(result, baseMod), modulus);
        }

        // Divide the exponent by 2
        exp = divideBigInts(exp, {true, {2}}).first;

        // Square the base and take modulo
        baseMod = modBigInts(multiplyBigInts(baseMod, baseMod), modulus);
    }

    return result;
}

// Đọc file test.inp
void readInput(const string& inputFile, int& x, int& y, BigInt& N, BigInt& e, vector<BigInt>& messages, vector<BigInt>& ciphertexts) {
    ifstream fin(inputFile);
    if (!fin.is_open()) {
        throw runtime_error("Cannot open file " + inputFile);
    }

    string line;
    getline(fin, line);
    stringstream ss(line);
    ss >> x >> y;

    if (x < 2 || x > 5){
        throw runtime_error("x must be in [2, 5]");
    }

    if (y < 10 || y > 20){
        throw runtime_error("y must be in [10, 20]");
    }

    // Đọc khóa công khai N và e
    getline(fin, line);
    size_t split = line.find(' ');
    N = hexToBigInt(line.substr(0, split));
    e = hexToBigInt(line.substr(split + 1));

    // Đọc tin nhắn
    for (int i = 0; i < x; ++i) {
        getline(fin, line);
        messages.push_back(hexToBigInt(line));
    }

    // Đọc bản mã
    for (int j = 0; j < y; ++j) {
        getline(fin, line);
        ciphertexts.push_back(hexToBigInt(line));
    }

    fin.close();
}

// Ghi file test.out
void writeOutput(const string& outputFile, const vector<int>& results) {
    ofstream fout(outputFile);
    if (!fout.is_open()) {
        throw runtime_error("Cannot open file " + outputFile);
    }

    for (size_t i = 0; i < results.size(); ++i) {
        if (i > 0) fout << " ";
        fout << results[i];
    }
    fout << endl;

    fout.close();
}

// Hàm chính
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    try {
        int x, y;
        BigInt N, e;
        vector<BigInt> messages, ciphertexts;

        auto start = std::chrono::high_resolution_clock::now();

        // Đọc dữ liệu
        readInput(inputFile, x, y, N, e, messages, ciphertexts);

        vector<int> results;

        // Xử lý từng tin nhắn
        for (const BigInt& message : messages) {
            BigInt cipher = modExp(message, e, N);
            int foundIndex = -1;

            // Tìm bản mã trong ciphertexts
            for (size_t j = 0; j < ciphertexts.size(); ++j) {
                if (compareBigInts(cipher, ciphertexts[j]) == 0) {
                    foundIndex = j;
                    break;
                }
            }

            results.push_back(foundIndex);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
        
        // Ghi kết quả ra file
        writeOutput(outputFile, results);

    } catch (const exception& ex) {
        cout << ex.what() << endl;
        return 1;
    }

    return 0;
}