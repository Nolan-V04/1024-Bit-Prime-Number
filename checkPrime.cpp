#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>

using BigInt = std::vector<uint32_t>;
const uint32_t BASE = 1000000000; // Base for storing digits (9 decimal digits)

// Converts a hexadecimal character to its integer value.
uint32_t hexToInt(char hexChar) {
    if (hexChar >= '0' && hexChar <= '9') return hexChar - '0';
    else if (hexChar >= 'A' && hexChar <= 'F') return hexChar - 'A' + 10;
    else if (hexChar >= 'a' && hexChar <= 'f') return hexChar - 'a' + 10;
    else throw std::invalid_argument("Invalid hex character");
}

// Multiplies a BigInt by 16.
void multiplyBy16(BigInt& bigInt) {
    uint64_t carry = 0;
    for (auto& digit : bigInt) {
        uint64_t product = static_cast<uint64_t>(digit) * 16 + carry;
        digit = product % BASE;
        carry = product / BASE;
    }
    if (carry > 0) bigInt.push_back(carry);
}

// Adds a single digit to a BigInt.
void addDigit(BigInt& bigInt, uint32_t digit) {
    uint64_t carry = digit;
    for (auto& element : bigInt) {
        uint64_t sum = static_cast<uint64_t>(element) + carry;
        element = sum % BASE;
        carry = sum / BASE;
        if (carry == 0) return;
    }
    if (carry > 0) bigInt.push_back(carry);
}

// Converts a hexadecimal string to a BigInt.
BigInt hexToBigInt(const std::string& hex) {
    if (hex.size() <= 8) { // Up to 8 hex digits can fit in uint32_t
        uint32_t smallNumber = stoul(hex, nullptr, 16);
        return BigInt{smallNumber};
    }
    BigInt result = {0};
    for (char hexChar : hex) {
        multiplyBy16(result);
        addDigit(result, hexToInt(hexChar));
    }
    return result;
}

// Compares two BigInts: returns 1 if a > b, -1 if a < b, and 0 if equal.
int compareBigInts(const BigInt& a, const BigInt& b) {
    if (a.size() > b.size()) return 1;
    if (a.size() < b.size()) return -1;
    for (size_t i = a.size(); i-- > 0;) {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

// Adds two BigInts.
BigInt addBigInts(const BigInt& a, const BigInt& b) {
    BigInt result;
    uint64_t carry = 0;
    size_t maxSize = std::max(a.size(), b.size());
    result.resize(maxSize);

    for (size_t i = 0; i < maxSize; ++i) {
        uint64_t sum = carry;
        if (i < a.size()) sum += a[i];
        if (i < b.size()) sum += b[i];
        result[i] = sum % BASE;
        carry = sum / BASE;
    }
    if (carry > 0) result.push_back(carry);

    return result;
}

// Subtracts one BigInt (b) from another (a).
BigInt subtractBigInts(const BigInt& a, const BigInt& b) {
    BigInt result = a;
    int64_t borrow = 0;

    for (size_t i = 0; i < b.size() || borrow > 0; ++i) {
        int64_t sub = borrow + (i < b.size() ? b[i] : 0);
        if (result[i] < sub) {
            result[i] += BASE - sub;
            borrow = 1;
        } else {
            result[i] -= sub;
            borrow = 0;
        }
    }

    while (result.size() > 1 && result.back() == 0) result.pop_back();

    return result;
}

// Multiplies two BigInts.
BigInt multiplyBigInts(const BigInt& a, const BigInt& b) {
    if (a.empty() || b.empty()) return {0};
    BigInt result(a.size() + b.size(), 0);

    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.size(); ++j) {
            uint64_t product = static_cast<uint64_t>(a[i]) * b[j] + result[i + j] + carry;
            result[i + j] = product % BASE;
            carry = product / BASE;
        }
        if (carry > 0) result[i + b.size()] += carry;
    }

    while (result.size() > 1 && result.back() == 0) result.pop_back();

    return result;
}

std::pair<BigInt, BigInt> divideBigInts(const BigInt& dividend, const BigInt& divisor) {
    if (compareBigInts(divisor, {0}) == 0) {
        throw std::invalid_argument("Division by zero");
    }

    if (compareBigInts(dividend, divisor) < 0) {
        return {{0}, dividend}; // Quotient = 0, Remainder = Dividend
    }

    BigInt quotient(dividend.size(), 0);
    BigInt remainder = dividend;

    size_t shift = dividend.size() - divisor.size();
    BigInt shiftedDivisor = divisor;
    shiftedDivisor.insert(shiftedDivisor.begin(), shift, 0);

    for (size_t i = shift + 1; i-- > 0;) {
        uint32_t q = 0;
        uint32_t low = 0, high = BASE - 1;

        // Binary search to find the maximum q such that (shiftedDivisor * q <= remainder)
        while (low <= high) {
            uint32_t mid = low + (high - low) / 2;
            BigInt testProduct = multiplyBigInts(shiftedDivisor, {mid});

            if (compareBigInts(testProduct, remainder) <= 0) {
                q = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }

        quotient[i] = q;
        remainder = subtractBigInts(remainder, multiplyBigInts(shiftedDivisor, {q}));

        if (i > 0) {
            shiftedDivisor.erase(shiftedDivisor.begin());
        }
    }

    while (quotient.size() > 1 && quotient.back() == 0) {
        quotient.pop_back();
    }
    while (remainder.size() > 1 && remainder.back() == 0) {
        remainder.pop_back();
    }

    return {quotient, remainder};
}

BigInt modMul(const BigInt& a, const BigInt& b, const BigInt& mod) {
    BigInt result = multiplyBigInts(a, b);
    return divideBigInts(result, mod).second;
}

BigInt modExp(BigInt base, BigInt exp, const BigInt& mod) {
    BigInt result = {1};
    base = divideBigInts(base, mod).second;

    while (compareBigInts(exp, {0}) > 0) {
        if (exp[0] % 2 != 0) {
            result = modMul(result, base, mod);
        }

        base = modMul(base, base, mod);

        BigInt temp;
        uint32_t carry = 0;
        for (int i = exp.size() - 1; i >= 0; --i) {
            uint32_t current = carry * BASE + exp[i];
            temp.insert(temp.begin(), current / 2);
            carry = current % 2;
        }
        exp = temp;

        while (exp.size() > 1 && exp.back() == 0) {
            exp.pop_back();
        }
    }

    return result;
}

std::random_device rd;
std::mt19937_64 gen(rd());

BigInt randomBigInt(const BigInt& n) {
    BigInt result(n.size(), 0);
    std::uniform_int_distribution<uint32_t> dist(0, BASE - 1);

    for (size_t i = 0; i < n.size(); ++i) {
        result[i] = dist(gen);
    }

    result[0] = std::max(result[0], uint32_t(2));
    if (compareBigInts(result, n) >= 0) {
        result = subtractBigInts(result, {2});
    }

    return result;
}

bool millerRabinTest(const BigInt& n, const BigInt& base) {
    BigInt nMinus1 = subtractBigInts(n, {1});
    BigInt d = nMinus1;
    size_t s = 0;

    while (d[0] % 2 == 0) {
        d = divideBigInts(d, {2}).first;
        ++s;
    }

    BigInt x = modExp(base, d, n);

    if (compareBigInts(x, {1}) == 0 || compareBigInts(x, nMinus1) == 0) {
        return true; // `x == 1` or `x == n-1`
    }

    for (size_t i = 0; i < s - 1; ++i) {
        x = modExp(x, {2}, n);
        if (compareBigInts(x, nMinus1) == 0) {
            return true;
        }
    }

    return false;
}

bool isPrime(const BigInt& n, size_t k = 5) {
    if (compareBigInts(n, {2}) < 0) return false;
    if (compareBigInts(n, {2}) == 0) return true;
    if (n[0] % 2 == 0) return false;

    for (size_t i = 0; i < k; ++i) {
        BigInt base = randomBigInt(n);
        if (!millerRabinTest(n, base)) {
            return false;
        }
    }

    return true;
}

// Read input from file
std::string readFile(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file " + filename + "! Retry.");
    }

    std::string line;
    fin >> line;
    fin.close();

    return line;
}

// Write result to file
void writeFile(const std::string& filename, int result) {
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Cannot open file " + filename + "! Retry.");
    }

    fout << result;
    fout.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    try {
    
        auto start = std::chrono::high_resolution_clock::now();

        // Read hexadecimal number from input file
        std::string hexNumber = readFile(inputFile);

        BigInt decNumber = hexToBigInt(hexNumber);
        int result = isPrime(decNumber) ? 1 : 0;

        // Write result to output file
        writeFile(outputFile, result);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}