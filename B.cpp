#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cassert>
#include <functional>
#include <complex>
#include <string>
#include <sstream>
#include <random>
#include <fstream>

const long double PI = std::acos(-1.0L);

struct UInt;
using ull = unsigned long long;
using ll = long long;
using ui = unsigned int;
using Number = UInt;

struct UInt { // Ukral otsuda https://github.com/dmkz/competitive-programming/blob/master/e-olymp.com/0317.cpp
    static const int BASE = (int)1e9; // ��������� ������� ���������
    static const int WIDTH = 9;       // ���������� ���������� ����, ������� �������� � ����� �����

    // ������ ��� ����� �����:
    std::vector<int> digits;

    // ������������
    UInt(int64_t number = 0);
    UInt(const std::string& s);
    UInt(const std::vector<int>& digits);

    // ������ ������������ � ���������:
    UInt& normalize(); // �������� ���������� ����� � �������� �� �������������� ���� ��������� [0, BASE)
    int compare(const UInt& other) const; // ��������� (������ = -1, ����� = 0, ������ = 1)

    // ������ ���������:
    UInt slow_mult(const UInt& other) const; // ��������� ������������ (�������� �������� ������ �� ������ ��������� �����)
    UInt fast_mult(const UInt& other) const; // ������� ������������ (�� ������ �������� �������������� ����� ����������� �����)
    UInt mult(const UInt& other) const; // ��������������� ����� ��������� �� ������ ����������������� ������

    // ����� �������:
    std::pair<UInt, UInt> div_mod(const UInt& other) const; // ����� ����� � ������� �� �������

    // ���������:
    UInt& operator+=(const int num);     // ����������� ���������
    UInt& operator+=(const UInt& other); // ����������� ��������
    UInt& operator-=(const int num);     // ��������� ���������
    UInt& operator-=(const UInt& other); // ��������� ��������
    UInt& operator*=(const int num);     // ��������� �� ��������
    UInt& operator*=(const UInt& other); // ��������� �� �������
    UInt& operator/=(const int num);     // ������� �� ��������
    UInt& operator/=(const UInt& other); // ������� �� �������
    UInt& operator%=(const UInt& other); // ������� �� ������� �� �������
};

std::istream& operator>>(std::istream&, UInt&); // ���� �� ������
std::ostream& operator<<(std::ostream&, const UInt&); // ����� � �����

UInt pow(UInt, int); // ���������� � �������
UInt gcd(UInt, UInt); // ���������� ����� ��������

UInt operator+(const UInt&, const UInt&);
UInt operator-(const UInt&, const UInt&);
UInt operator*(const UInt&, const UInt&);
UInt operator/(const UInt&, const UInt&);
UInt operator%(const UInt&, const UInt&);

UInt operator+(const UInt&, const int);
UInt operator+(const int, const UInt&);
UInt operator-(const UInt&, const int);
UInt operator*(const UInt&, const int);
UInt operator*(const int, const UInt&);
UInt operator/(const UInt&, const int);
UInt operator^(const UInt&, const int); // ���������� � �������

bool operator<(const UInt&, const UInt&);
bool operator>(const UInt&, const UInt&);
bool operator<=(const UInt&, const UInt&);
bool operator>=(const UInt&, const UInt&);
bool operator==(const UInt&, const UInt&);
bool operator!=(const UInt&, const UInt&);

UInt& UInt::normalize() {
    if (digits.size() == 0)
        throw new std::bad_exception;

    while (digits.back() == 0 && (int)digits.size() > 1) digits.pop_back();
    for (auto d : digits) assert(0 <= d && d < BASE);
    return *this;
}

// ����������� �� ��������� ������
UInt::UInt(int64_t number) {
    assert(number >= 0);
    do {
        digits.push_back(number % BASE);
        number /= BASE;
    } while (number > 0);
    normalize();
}

// ����������� �� ������� �� ����:
UInt::UInt(const std::vector<int>& digits) : digits(digits) {
    normalize();
}

// ����������� �� �������:
UInt::UInt(const std::string& s) {
    const int size = (int)s.size();
    for (int idGroup = 1, nGroups = size / WIDTH; idGroup <= nGroups; ++idGroup) {
        digits.push_back(std::stoi(s.substr(size - idGroup * WIDTH, WIDTH)));
    }
    if (size % WIDTH != 0) {
        digits.push_back(std::stoi(s.substr(0, size % WIDTH)));
    }
    normalize();
}

// ����������� ���������:
UInt& UInt::operator+=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this += UInt(num);
    }
    int rem = num;
    for (int i = 0; rem > 0; ++i) {
        if (i >= (int)digits.size()) digits.push_back(0);
        rem += digits[i];
        if (rem >= BASE) {
            digits[i] = rem - BASE;
            rem = 1;
        }
        else {
            digits[i] = rem;
            rem = 0;
        }
    }
    return this->normalize();
}

// ����������� ��������:
UInt& UInt::operator+=(const UInt& other) {
    if (other.digits.size() == 1u) {
        return *this += other.digits[0];
    }
    const int s1 = this->digits.size();
    const int s2 = other.digits.size();
    int rem = 0;
    for (int i = 0; i < s1 || i < s2 || rem > 0; ++i) {
        int d1 = i < s1 ? this->digits[i] : (digits.push_back(0), 0);
        int d2 = i < s2 ? other.digits[i] : 0;
        rem += d1 + d2;
        auto div = rem / BASE;
        digits[i] = rem - div * BASE;
        rem = div;
    }
    return this->normalize();
}

// ��������� ���������:
UInt& UInt::operator-=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this -= UInt(num);
    }
    int rem = -num;
    for (int i = 0; i < (int)digits.size() && rem < 0; ++i) {
        rem += digits[i];
        if (rem < 0) { // �������� ������
            digits[i] = rem + BASE;
            rem = -1;
        }
        else {
            digits[i] = rem;
            rem = 0;
        }
    }
    assert(rem == 0);
    return this->normalize();
}

// ��������� ��������:
UInt& UInt::operator-=(const UInt& other) {
    if (other.digits.size() == 1u) {
        return *this -= other.digits[0];
    }
    const int s1 = this->digits.size();
    const int s2 = other.digits.size();
    assert(s1 >= s2);
    int rem = 0;
    for (int i = 0; i < s1; ++i) {
        int d2 = i < s2 ? other.digits[i] : 0;
        rem += this->digits[i] - d2;
        if (rem < 0) {
            digits[i] = rem + BASE;
            rem = -1;
        }
        else {
            digits[i] = rem;
            rem = 0;
            if (i >= s2) break;
        }
    }
    assert(rem == 0); // ����� *this < other
    return this->normalize();
}

// ��������� �� ��������:
UInt& UInt::operator*=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this *= UInt(num);
    }
    int64_t rem = 0;
    for (auto& d : digits) {
        rem += 1LL * d * num;
        auto div = rem / BASE;
        d = rem - div * BASE;
        rem = div;
    }
    if (rem > 0) digits.push_back(rem);
    return this->normalize();
}

// ��������� ������������:
UInt UInt::slow_mult(const UInt& other) const {
    if (other.digits.size() == 1u) {
        return *this * other.digits[0];
    }
    const int s1 = (int)this->digits.size();
    const int s2 = (int)other.digits.size();
    std::vector<int> temp(s1 + s2);
    for (int i = 0; i < s1; ++i) {
        int64_t rem = 0;
        for (int j = 0; j < s2; ++j) {
            rem += temp[i + j] + 1LL * this->digits[i] * other.digits[j];
            auto div = rem / BASE;
            temp[i + j] = rem - div * BASE;
            rem = div;
        }
        if (rem > 0) temp[i + s2] += rem;
    }
    return UInt(temp);
}

// ������� ��������� �� ������ �������� �������������� �����:
UInt UInt::fast_mult(const UInt& other) const {
    if (other.digits.size() == 1u) {
        return *this * other.digits[0];
    }

    // �������� ����� � ����� num:
    std::function<int(int, int)> reverse = [](int number, int nBits) {
        int res = 0;
        for (int i = 0; i < nBits; ++i) {
            if (number & (1 << i)) {
                res |= 1 << (nBits - 1 - i);
            }
        }
        return res;
    };

    typedef std::complex<long double> complex;
    // ������� �������������� �����:
    std::function<void(std::vector<complex>&, bool)> fft = [&reverse](std::vector<complex>& a, bool invert) {
        const int n = (int)a.size();
        int nBits = 0;
        while ((1 << nBits) < n) ++nBits;

        for (int i = 0; i < n; ++i) {
            if (i < reverse(i, nBits)) {
                std::swap(a[i], a[reverse(i, nBits)]);
            }
        }

        for (int len = 2; len <= n; len <<= 1) {
            auto ang = 2 * PI / len * (invert ? -1 : 1);
            complex wlen(std::cos(ang), std::sin(ang));
            for (int i = 0; i < n; i += len) {
                complex w(1);
                for (int j = 0; j < len / 2; ++j) {
                    complex u = a[i + j];
                    complex v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) {
            for (int i = 0; i < n; ++i) {
                a[i] /= n;
            }
        }
    };

    // �������������� ������� �� ����������� ������������� fa � fb:
    // ��� ��� ���������� ������ �������� ��-�� ���������� � ��������� ������, ��������� ������� ���������� ��������:
    assert(BASE == 1000 * 1000 * 1000);
    std::function<std::vector<complex>(const UInt&)> prepare = [](const UInt& number) {
        std::vector<complex> result;
        result.reserve(3 * number.digits.size());
        for (auto d : number.digits) {
            result.push_back(d % 1000);
            result.push_back(d / 1000 % 1000);
            result.push_back(d / 1000000);
        }
        return result;
    };

    auto fa = prepare(*this);
    auto fb = prepare(other);

    // ��������� ������ �������� �� ��������� ������� ������:
    int n = 1;
    while (n < (int)std::max(fa.size(), fb.size())) n *= 2;
    n *= 2;
    fa.resize(n);
    fb.resize(n);

    // �������� ������ �������������� �����:
    fft(fa, false);
    fft(fb, false);
    // ����������� ����������:
    for (int i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }
    // �������� �������� �������������� �����:
    fft(fa, true);
    // �������� ����� � ������������:
    std::vector<int64_t> temp(n);
    for (int i = 0; i < (int)fa.size(); ++i) {
        temp[i] = int64_t(fa[i].real() + 0.5);
    }
    // �� �������� ��� �������� � ������� �������:
    int64_t carry = 0;
    for (int i = 0; i < n || carry > 0; ++i) {
        if (i >= n) temp.push_back(0);
        temp[i] += carry;
        carry = temp[i] / 1000;
        temp[i] -= carry * 1000;
        assert(temp[i] >= 0);
    }
    // ��������� �����:
    std::vector<int> res;
    res.reserve(this->digits.size() + other.digits.size());

    for (int i = 0; i < n; i += 3) {
        int c = temp[i];
        int b = i + 1 < n ? temp[i + 1] : 0;
        int a = i + 2 < n ? temp[i + 2] : 0;
        res.push_back(c + 1000 * (b + 1000 * a));
    }
    return UInt(res);
}

// ��������������� ����� ���������:
UInt UInt::mult(const UInt& other) const {
    // ����� ������ ���������:
    int len1 = (int)this->digits.size();
    int len2 = (int)other.digits.size();
    int temp = 3 * std::max(len1, len2);
    int pow = 1;
    while (pow < temp) pow *= 2;
    pow *= 2;
    int op1 = len1 * len2;
    int op2 = 3 * pow * std::log(pow) / std::log(2);
    return op1 >= 15 * op2 ? fast_mult(other) : slow_mult(other);
}

// ������� �� ��������:
UInt& UInt::operator/=(const int num) {
    assert(num > 0);
    if (num >= BASE) {
        return *this /= UInt(num);
    }
    int64_t rem = 0;
    for (int j = (int)digits.size() - 1; j >= 0; --j) {
        (rem *= BASE) += digits[j];
        auto div = rem / num;
        digits[j] = div;
        rem -= div * num;
    }
    return this->normalize();
}

// ������� �� ������� �� ��������:
int operator%(const UInt& a, const int num) {
    assert(num > 0);
    int64_t rem = 0;
    for (int i = (int)a.digits.size() - 1; i >= 0; --i) {
        ((rem *= UInt::BASE) += a.digits[i]) %= num;
    }
    return rem;
}

// ����� ����� � ������� �� �������:
std::pair<UInt, UInt> UInt::div_mod(const UInt& other) const {
    if (other.digits.size() == 1u) {
        return { std::move(*this / other.digits[0]), *this % other.digits[0] };
    }
    const int norm = BASE / (other.digits.back() + 1);
    const UInt a = *this * norm;
    const UInt b = other * norm;
    const int a_size = (int)a.digits.size();
    const int b_size = (int)b.digits.size();
    UInt q, r;
    q.digits.resize(a_size);
    for (int i = a_size - 1; i >= 0; --i) {
        r *= BASE;
        r += a.digits[i];
        int s1 = (int)r.digits.size() <= b_size ? 0 : r.digits[b_size];
        int s2 = (int)r.digits.size() <= b_size - 1 ? 0 : r.digits[b_size - 1];
        int d = (1LL * BASE * s1 + s2) / b.digits.back();
        auto temp = b * d;
        while (r < temp) {
            r += b;
            --d;
        }
        r -= temp;
        q.digits[i] = d;
    }
    return { std::move(q.normalize()), std::move(r /= norm) };
}

// ���������: result < 0 (������), result == 0 (�����), result > 0 (������)
int UInt::compare(const UInt& other) const {
    if (this->digits.size() > other.digits.size()) return 1;
    if (this->digits.size() < other.digits.size()) return -1;
    for (int i = (int)digits.size() - 1; i >= 0; --i) {
        if (this->digits[i] > other.digits[i]) return 1;
        if (this->digits[i] < other.digits[i]) return -1;
    }
    return 0;
}

// ��������� ���������:
bool operator< (const UInt& a, const UInt& b) { return a.compare(b) < 0; }
bool operator> (const UInt& a, const UInt& b) { return a.compare(b) > 0; }
bool operator==(const UInt& a, const UInt& b) { return a.compare(b) == 0; }
bool operator<=(const UInt& a, const UInt& b) { return a.compare(b) <= 0; }
bool operator>=(const UInt& a, const UInt& b) { return a.compare(b) >= 0; }
bool operator!=(const UInt& a, const UInt& b) { return a.compare(b) != 0; }

// ���� �� ������:
std::istream& operator>>(std::istream& is, UInt& number) {
    std::string s;
    is >> s;
    number = UInt(s);
    return is;
}

// ����� � �����:
std::ostream& operator<<(std::ostream& os, const UInt& number) {
    os << number.digits.back();
    for (int i = (int)number.digits.size() - 2; i >= 0; --i) {
        os << std::setw(UInt::WIDTH) << std::setfill('0') << number.digits[i];
    }
    return os << std::setfill(' ');
}

// �����:
UInt operator+(const UInt& a, const UInt& b) {
    return UInt(a) += b;
}

// ��������:
UInt operator-(const UInt& a, const UInt& b) {
    return UInt(a) -= b;
}

// ������������:
UInt operator*(const UInt& a, const UInt& b) {
    return a.mult(b);
}

// �������:
UInt operator/(const UInt& a, const UInt& b) {
    return a.div_mod(b).first;
}

// ������ �������:
UInt operator%(const UInt& a, const UInt& b) {
    return a.div_mod(b).second;
}

// ���������:
UInt& UInt::operator*=(const UInt& other) {
    return *this = *this * other;
}

// ������� � �������������:
UInt& UInt::operator/=(const UInt& other) {
    return *this = *this / other;
}

// ������ ������� � �������������:
UInt& UInt::operator%=(const UInt& other) {
    return *this = *this % other;
}

UInt operator+(const UInt& a, const int b) { return UInt(a) += b; }
UInt operator+(const int a, const UInt& b) { return b * a; }
UInt operator-(const UInt& a, const int b) { return UInt(a) -= b; }
UInt operator*(const UInt& a, const int b) { return UInt(a) *= b; }
UInt operator*(const int a, const UInt& b) { return b * a; }
UInt operator/(const UInt& a, const int b) { return UInt(a) /= b; }
UInt operator^(const UInt& a, const int n) { return pow(a, n); } // ���������� � �������

// ���������� � �������:
UInt pow(UInt a, int n) {
    UInt res = 1;
    while (n > 0) {
        if (n % 2 != 0) res *= a;
        a *= a;
        n /= 2;
    }
    return res;
}

// ���������� ����� ��������:
UInt gcd(UInt a, UInt b) {
    while (b != 0) {
        auto rem = a % b;
        a = b;
        b = rem;
    }
    return a;
}


Number get_random(Number a, Number b)
{
    static std::mt19937 mt(std::time(0));

    return mt() % (b - a) + a;
}

class GroupZp
{
public:
    GroupZp()
        : value(1) {}

    GroupZp(Number val)
        : value(val) {}

    GroupZp operator*(const GroupZp& r)
    {
        return GroupZp((value * r.value) % order);
    }

    GroupZp operator/(const GroupZp& r)
    {
        return GroupZp((value * r.inverse().value) % order);
    }

    GroupZp inverse() const
    {
        return GroupZp(power(order-2));
    }

    GroupZp power(Number p) const
    {
        if (p == 1)
            return GroupZp(Number(value));
        if (p == 0)
            return GroupZp(Number(1));

        GroupZp power_half = power(p / 2);

        if (p % 2)
            return power_half * power_half * GroupZp(value);
        else
            return power_half * power_half;
    }

    static void set_order(Number ord)
    {
        order = ord;
    }

    Number get_value() const
    {
        return value;
    }

    void set_value(Number val)
    {
        value = val;
    }

private:
    Number value;
    static Number order;
};

Number GroupZp::order = 2;

std::istream& operator>>(std::istream& in, GroupZp& r)
{
    Number n;
    in >> n;

    r.set_value(n);
    return in;
}

std::ostream& operator<<(std::ostream& out, const GroupZp& r)
{
    out << r.get_value();
    return out;
}

Number char_to_number(char symbol) 
{
    if (symbol >= 48 && symbol <= 57)
        return symbol - 48;
    if (symbol >= 65 && symbol <= 90)
        return symbol - 55;
    if (symbol >= 97 && symbol <= 122)
        return symbol - 61;
    if (symbol == 32)
        return 62;
    if (symbol == 46)
        return 63;
    return 64;
}
char number_to_char(int num)
{
    if (num >= 0 && num <= 9)
        return (num + 48);
    if (num >= 10 && num <= 35)
        return (num + 55);
    if (num >= 36 && num <= 61)
        return (num + 61);
    if (num == 62)
        return (32);
    if (num == 63)
        return (46);
    return 64;
}

Number char_to_number_fixed(char symbol) 
{
    return char_to_number(symbol);
}
char number_to_char_fixed(int num)
{
    return number_to_char(num);
}

GroupZp elgamal_encode(Number c, GroupZp g, GroupZp k, Number b)
{
    return k.power(b) * GroupZp(c);
}

Number elgamal_decode(GroupZp n, Number a, GroupZp gb)
{
    GroupZp k = gb.power(a);
    GroupZp h = n * k.inverse();
    auto r = h.get_value();
    return  h.get_value();
}

std::vector<Number> change_base(const std::vector<Number>& num, Number from, Number to)
{
    Number a = 0;
    Number p = 1;

    for (int i = 0; i < num.size(); ++i)
    {
        a += p * num[i];
        p *= from;
    }

    std::vector<Number> res;

    while (a > 0)
    {
        res.push_back(a % to);
        a /= to;
    }
    return res;
}

int main_encode()
{
    Number p;
    std::cin >> p;

    GroupZp::set_order(p);
    GroupZp g, k;
    std::cin >> g >> k;


    std::string s;
    std::getline(std::cin, s);
    std::getline(std::cin, s);
    std::vector<Number> num;
    for (int i = 0; i < s.size(); ++i)
    {
        num.push_back(char_to_number_fixed(s[i]));
    }

    auto num1 = change_base(num, 64, p);

    for (Number n : num1)
    {
        Number b = get_random(1, p);
        GroupZp elem = elgamal_encode(n, g, k, b);
        std::cout << g.power(b) << ' ' << elem << std::endl;
    }
    return 0;
}

int main_decode()
{
    Number p;
    std::cin >> p;

    GroupZp::set_order(p);

    Number k;
    std::cin >> k;

    std::vector<Number> num;

    while (std::cin)
    {
        GroupZp gb, h;
        try
        {
            std::cin >> gb >> h;
        }
        catch (...)
        {
            break;
        }

        Number c = elgamal_decode(h, k, gb);
        num.push_back(c);
    }

    std::vector<Number> num1 = change_base(num, p, 64);

    for (int i = 0; i < num1.size(); ++i)
    {
        std::stringstream converter;
        converter << num1[i];
        int converted;
        converter >> converted;
        std::cout << number_to_char(converted);
    }

    return 0;
}

int main()
{
    if (0)
    {
        return main_encode();
    }
    else
    {
        return main_decode();
    }
}