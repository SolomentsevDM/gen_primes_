#include <iostream>
#include<vector>
#include <fstream>
#include<clocale>
#include <thread>
#include<mutex>
#include<complex>
#include<ctime>
#include<chrono>
typedef std::complex<double> base;
# define PI           3.14159265358979323846

std::vector<int>mod_pow(std::vector<int> n, std::vector<int> deg, std::vector<int> m);
std::vector<int> vec_pow(std::vector<int> num, std::vector<int> dig);
std::vector<int> mod(std::vector<int> num1, std::vector<int>num2);
std::vector<int> operator / (std::vector<int> num1, std::vector<int>num2);
std::vector<int> operator *(std::vector<int> num1, std::vector<int> num2);
std::vector<int> operator *(std::vector<int> num, int dig);
std::vector<int> operator + (std::vector<int> num1, int dig);
std::vector<int> operator + (std::vector<int> num1, std::vector<int>num2);
std::vector<int> operator - (std::vector<int> num1, std::vector<int>num2);
bool operator >=(std::vector<int> num1, std::vector<int>num2);
bool operator <=(std::vector<int> num1, std::vector<int>num2);
bool operator >(std::vector<int> num1, std::vector<int>num2);
bool operator <(std::vector<int> num1, std::vector<int>num2);
std::string print_vec(std::vector<int> vec);
bool test_M_R(std::vector<int> n, int thread_num);
void get_prime(std::vector<int> start, std::vector<int> end, int thread_num);
void printFile(std::vector<int> vec);
std::vector<int> get_rnd(int size, int thread_num, std::vector<int> max);
void normalize(std::vector<int>& l);
void fft(std::vector<base>& a, bool invert);
void multiply(const std::vector<int>& a, const std::vector<int>& b, std::vector<int>& res);
void progress_bar();
std::mutex mtx;
int k = 0;
int in_K = -1;
bool primes_find = false;
int thread_number = 4;
int thread_complete = 0;
std::vector<int> thread_n = { 4 };
std::thread* thread_arr = new std::thread[thread_number];
std::ofstream file;
int main()
{
    setlocale(LC_ALL, "Russian");
    while (true)
    {
        std::cout << "Введите количество простых чисел: " << std::flush;
        if ((std::cin >> in_K).good() && in_K > 0) break;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cout << "Неверный ввод, повторите.\n";
        }
        else
        {
            std::cout << "Количество должно быть больше нуля\n";
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    int deg_10 = -1;
    while (true)
    {
        std::cout << "Порядок числа.\n" << std::flush;
        std::cout << "Пример: \nN\t\t\tИнтервал генерации простых чисел\n3\t\t\t(1000, 10000)\n6\t\t\t(1000000, 10000000)\n9\t\t\t(1000000000, 10000000000)\n";
        std::cout << "Введите порядок числа: " << std::flush;
        if ((std::cin >> deg_10).good() && deg_10 > 0 && deg_10 < 80) break;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cout << "Неверный ввод, повторите.\n";
        }
        else
        {
            std::cout << "Значение порчдка должно быть натуральным числом и меньше 80\n";
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    time_t sec = time(NULL);
    tm* timeinfo = localtime(&sec);
  
    std::string new_name = "primes--" + std::to_string(2000+timeinfo->tm_year - 100) +"-" + std::to_string(timeinfo->tm_mon) +"-"+ std::to_string(timeinfo->tm_mday) + "--"
        + std::to_string(timeinfo->tm_hour) +"-"+ std::to_string(timeinfo->tm_min)+"-" + std::to_string(timeinfo->tm_sec) + +".txt";
    file.open(new_name);
    std::vector<int> num;
    std::thread progress;
    for (int i = 0; i < deg_10 + 1; i++)
    {
        num.push_back(0);
    }
    num[num.size() - 1] = 1;
    std::vector<int> step = (num * 10 - num) / thread_n;
    for (int i = 0; i < thread_number; i++)
    {
        thread_arr[i] = (std::thread(get_prime, num + step * i, num + step * (i + 1), i));
    }
    progress = (std::thread(progress_bar));
    for (int i = 0; i < thread_number; i++)
    {
        thread_arr[i].join();
    }
    progress.join();
    file.close();
    return 0;
}
void progress_bar()
{
    int start = 0;
    while (thread_complete!=thread_number)
    {
        if (start != k)
        {
            for (int i = 0; i < 20 / in_K; i++)
                std::cout << "#";
            start++;
        }
        for (int i = 0; i < 20; i++)
        {
            std::cout << "#";
            _sleep(300);
        }
        _sleep(200);
        for (int i = 0; i < 20; i++)
        {
            printf("\b \b");
        }
        _sleep(400);
    }
    std::cout << "\nComplete!";
}
void printFile(std::vector<int> vec)
{
    mtx.lock();
    file << print_vec(vec) << "\n\n";
    k++;
    mtx.unlock();
}
void get_prime(std::vector<int> start, std::vector<int> end, int thread_num)
{
    std::vector<int> one = { 1 };
    end = end - one;
    start = get_rnd(end.size(), thread_num,end);
    if (start[0] % 2 == 0) start[0] = start[0] + 1;
    while (true)
    {
        if (primes_find) break;
        for (; start < end; start = start + 2)
        {
            if (primes_find) {
                thread_complete++;
                break;
            }
            if (start[0] != 5)
            {
                if (test_M_R(start, thread_num))
                {
                    if (k < in_K)
                    {
                        printFile(start);
                        start = get_rnd(end.size(), thread_num,end);
                        if (start[0] % 2 == 0) start[0] = start[0] + 1;
                        if (k == in_K)
                        {
                            primes_find = true;
                        }
                    }
                    else
                    {
                        primes_find = true;
                        thread_complete++;
                        break;
                    }
                }
            }
        }
        start = get_rnd(end.size(), thread_num, end);
        if (start[0] % 2 == 0) start[0] = start[0] + 1;

    }
}
std::vector<int> get_rnd(int size, int thread_num, std::vector<int> max)
{
    mtx.lock();
    srand(time(0) * (thread_num + 1));
    std::vector<int> vec(size);
    do
    {

        for (int i = 0; i < vec.size(); i++)
        {
            if (i == vec.size() - 1)
            {
                vec[i] = rand() % max[i] + 1;
            }
            else
            {
                vec[i] = rand() % 10;
            }
        }
    } while (vec > max);
    mtx.unlock();
    return vec;
}
bool pre_check(std::vector<int> n)
{
    std::vector<int> dig = { 7,5,2 };
    if (n[0] % 2 == 0) std::cout << "hahahah";
    for (std::vector<int> i = { 3 }; i < dig; i = i + 2)
    {
        std::vector<int> modd = mod(n, i);
        if (modd.size() == 1 && modd[0] == 0) return false;
    }
    return true;
}
bool test_M_R(std::vector<int> n, int thread_num)
{
    if (primes_find) return false;

    if (!pre_check(n)) return false;

    std::vector<int> one = { 1 };
    std::vector<int> two = { 2 };
    std::vector<int> s = { 0 };
    std::vector<int> modd = mod(n - one, { 2 });
    std::vector<int> d = n - one;
    while (modd.size() == 1 && modd[0] == 0)
    {
        if (primes_find) break;
        modd = mod(d, { 2 });
        if (modd.size() == 1 && modd[0] == 0)
        {
            d = d / two;
            s = s + 1;
        }
    }
    for (int i = 0; i < 5; i++)
    {
        if (primes_find) break;
        int size = 1 + rand() % (n.size());
        if (size != 1) size--;
        std::vector<int> A(size);
        srand(time(0) * (thread_num + 1));
        mtx.lock();
        do
        {
            for (int i = 0; i < A.size(); i++)
            {
                if (i == A.size() - 1)
                {
                    A[i] = rand() % 9 + 1;
                }
                else
                {
                    A[i] = rand() % 10;
                }
            }
        } while (A >= n);
        mtx.unlock();
        std::vector<int> res = mod_pow(A, d, n);
        if ((res.size() == 1 && res[0] == 1) || res == n - one) continue;
        for (std::vector<int> r = { 1 }; r < s; r = r + 1)
        {
            if (primes_find) break;
            res = mod_pow(res, { 2 }, n);
            if (res == one) return false;
            if (res == n - one) break;
        }
        if (!(res == n - one)) return false;
    }

    return true;
}
std::string print_vec(std::vector<int> vec)
{
    std::string str;
    for (int i = vec.size() - 1; i >= 0; i--)
    {
        str += vec[i] + '0';
    }
    return str;
}
std::vector<int>mod_pow(std::vector<int> n, std::vector<int> deg, std::vector<int> m)
{
    if (primes_find) return { 1 };
    if (deg.size() == 1 && deg[0] == 0) return { 1 };
    std::vector<int> two = { 2 };
    std::vector<int> z = mod_pow(n, deg / two, m);
    if (deg[0] % 2 == 0)
    {
        multiply(z, z, z);
        return mod(z, m);
    }
    else
    {
        multiply(z, z, z);
        multiply(z, n, z);
        return mod(z, m);
    }

}
//std::vector<int> vec_pow(std::vector<int> num, std::vector<int> dig)
//{
//    std::vector<int> res = { 1 };
//    std::vector<int> n = { 2 };
//    while (!(dig.size() == 1 && dig[0] == 0))
//    {
//        if (mod(dig, n)[0] == 1)
//        {
//            res = res * num;
//        }
//        num = num * num;
//        dig = dig / n;
//    }
//    return res;
//}
std::vector<int> operator / (std::vector<int> num1, std::vector<int>num2)
{
    /*std::vector<int> Q={0};
    std::vector<int> R = num1;
    while (R > num2)
    {
        R = R - num2;
        Q = Q + 1;
    }
    return Q;*/
    if (num1 < num2) return std::vector<int>{0};
    std::vector<int> tmp;
    int k = 0;
    for (int i = num1.size() - 1; k < num2.size(); i--, k++)
    {
        tmp.insert(tmp.begin(), num1[i]);
    }
    int ind = num1.size() - 1 - k;
    while (tmp < num2)
    {
        tmp.insert(tmp.begin(), num1[ind]);
        ind--;
    }
    num1.erase(num1.end() - tmp.size(), num1.end());
    std::vector<int>del = num2;
    std::vector<int> quotient;
    int times = 1;
    while (del <= tmp)
    {
        del = num2;
        times = times + 1;
        del = del * times;
    }
    quotient.push_back(times - 1);
    tmp = tmp - quotient * num2;
    while (num1.size() != 0)
    {
        tmp.insert(tmp.begin(), num1[num1.size() - 1]);
        for (int i = tmp.size() - 1; i >= 1 && tmp[i] == 0; i--)
        {
            tmp.erase(tmp.begin() + i);
        }
        num1.erase(num1.end() - 1);
        if (tmp >= num2)
        {
            times = 1;
            del = num2;
            while (del <= tmp)
            {
                del = num2;
                times = times + 1;
                del = del * times;
            }
            quotient.insert(quotient.begin(), times - 1);
            tmp = tmp - num2 * (times - 1);
        }
        else
        {
            quotient.insert(quotient.begin(), 0);
        }
    }
    return quotient;
}
std::vector<int> mod(std::vector<int> num1, std::vector<int>num2)
{
    if (num1 < num2) return num1;
    std::vector<int> t = num1 - (num1 / num2) * num2;
    return t;
}
std::vector<int> operator *(std::vector<int> num1, std::vector<int> num2)
{
    if (num2.size() == 1 && num2[0] == 1) return num1;
    if (num1.size() == 1 && num1[0] == 1) return num2;
    std::vector<int> res;
    bool flag = false;
    int size = std::max(num1.size(), num2.size());
    int times;
    int p;
    for (int i = 0; i < num2.size(); i++)
    {
        std::vector<int> tmp_res;
        flag = false;
        for (int j = 0; j < num1.size(); j++)
        {
            times = num2[num2.size() - i - 1] * num1[j];
            if (flag) times += p;
            if (times >= 10)
            {
                flag = true;
                p = times / 10;
            }
            else flag = false;
            tmp_res.push_back(times % 10);

        }
        if (flag) tmp_res.push_back(p);
        if (res.empty())
        {
            res = tmp_res;
        }
        else
        {
            res.insert(res.begin(), 0);
            res = res + tmp_res;
        }
    }
    if (flag)
    {

    }
    return res;
}
std::vector<int> operator *(std::vector<int> num, int dig)
{
    std::vector<int> res;
    bool flag = false;
    int size = num.size();
    int times;
    int p;
    for (int i = 0; i < size; i++)
    {
        std::vector<int> tmp_res;
        flag = false;
        times = num[size - i - 1] * dig;
        if (flag) times += p;
        if (times >= 10)
        {
            flag = true;
            p = times / 10;
        }
        else flag = false;
        tmp_res.push_back(times % 10);
        if (flag) tmp_res.push_back(p);
        if (res.empty())
        {
            res = tmp_res;
        }
        else
        {
            res.insert(res.begin(), 0);
            res = res + tmp_res;
        }
    }
    if (flag)
    {

    }
    for (int i = res.size() - 1; i > 0 && res[i] == 0; i--)
    {
        res.erase(res.begin());
    }
    //std::reverse(res.num.begin(), res.num.end());
    return res;
}
std::vector<int> operator + (std::vector<int> num1, int dig)
{
    num1[0] += dig;
    std::vector<int> res;
    bool flag = false;
    int t = 0;
    for (int i = 0; i < num1.size(); i++)
    {
        if (flag)
        {
            num1[i] += t;
            flag = false;
        }
        if (num1[i] >= 10)
        {
            res.push_back(num1[i] % 10);
            flag = true;
            t = num1[i] / 10;
        }
        if (!flag)
        {
            res.push_back(num1[i]);
        }
    }
    if (flag)
    {
        res.push_back(t);
    }
    if (res.size() == 0) return num1;
    return res;
}
std::vector<int> operator + (std::vector<int> num1, std::vector<int>num2)
{
    std::vector<int> res;
    bool flag = false;
    int size = std::max(num1.size(), num2.size());
    int p;
    int min = std::min(num1.size(), num2.size());
    int plus;
    if (num1.size() < num2.size())
    {
        std::swap(num1, num2);
    }
    for (int i = 0; i < size; i++)
    {
        if (i >= min)
        {
            plus = num1[i];
        }
        else
        {
            plus = (num1[i]) + num2[i];
        }
        if (flag) plus += p;
        if (plus >= 10)
        {
            flag = true;
            p = plus / 10;
        }
        else flag = false;
        res.push_back(plus % 10);
    }
    if (flag)
    {
        res.push_back(1);
    }
    return res;
}
std::vector<int> operator - (std::vector<int> num1, std::vector<int>num2)
{
    int size = std::max(num1.size(), num2.size());
    int min = std::min(num1.size(), num2.size());
    std::vector<int> res(size);
    int plus;
    for (int i = 0; i < size; i++)
    {
        if (i >= min)
        {
            res[i] = num1[i];
        }
        else
        {
            res[i] = (num1[i]) - num2[i];
        }
    }
    for (int i = 0; i < res.size() - 1; i++)
    {
        if (res[i] < 0) { //если меньше - заем
            int carryover = (res[i] + 1) / 10 - 1;
            res[i + 1] += carryover;
            res[i] -= carryover * 10;
        }
    }
    for (int i = res.size() - 1; i >= 1 && res[i] == 0; i--)
    {
        res.erase(res.begin() + i);
    }
    return res;
}
bool operator >=(std::vector<int> num1, std::vector<int>num2)
{
    if (num1.size() > num2.size())return true;
    if (num1.size() < num2.size())return false;
    for (int i = num1.size() - 1; i >= 0; i--)
    {
        if (num1[i] > num2[i]) return true;
        if (num1[i] < num2[i]) return false;

    }
    return true;
}
bool operator <=(std::vector<int> num1, std::vector<int>num2)
{
    if (num1.size() > num2.size())return false;
    if (num1.size() < num2.size())return true;
    for (int i = num1.size() - 1; i >= 0; i--)
    {
        if (num1[i] < num2[i]) return true;
        if (num1[i] > num2[i]) return false;
    }
    return true;
}
bool operator >(std::vector<int> num1, std::vector<int>num2)
{
    if (num1.size() > num2.size())return true;
    if (num1.size() < num2.size())return false;
    for (int i = num1.size() - 1; i >= 0; i--)
    {
        if (num1[i] > num2[i]) return true;
        if (num1[i] < num2[i]) return false;
    }
    return false;
}
bool operator <(std::vector<int> num1, std::vector<int>num2)
{
    if (num1.size() > num2.size())return false;
    if (num1.size() < num2.size())return true;
    for (int i = num1.size() - 1; i >= 0; i--)
    {
        if (num1[i] < num2[i]) return true;
        if (num1[i] > num2[i]) return false;
    }
    return false;
}
bool operator ==(std::vector<int> num1, std::vector<int>num2)
{
    if (num1.size() != num2.size()) return false;
    for (int i = 0; i < num1.size(); i++)
    {
        if (num1[i] != num2[i]) return false;
    }
    return true;
}
std::vector<int> vec_pow(std::vector<int> num, std::vector<int> dig)
{
    std::vector<int> res = { 1 };
    std::vector<int> n = { 2 };
    while (!(dig.size() == 1 && dig[0] == 0))
    {
        if (mod(dig, n)[0] == 1)
        {
            multiply(res, num, res);
        }
        multiply(num, num, num);
        dig = dig / n;
    }
    return res;
}
void multiply(const std::vector<int>& a, const std::vector<int>& b, std::vector<int>& res) {
    int max = std::max(a.size(), b.size());
    if (max <= 10) { //если число короче то применять наивное умножение
        std::vector<int>result(2 * max);
        for (int i = 0; i < a.size(); ++i)
            for (int j = 0; j < b.size(); ++j) {
                result[i + j] += a[i] * b[j];
            }
        res = result;
    }
    else
    {
        std::vector<base> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        size_t n = 1;
        while (n < std::max(a.size(), b.size()))  n <<= 1;
        n <<= 1;
        fa.resize(n), fb.resize(n);

        fft(fa, false), fft(fb, false);
        for (size_t i = 0; i < n; ++i)
            fa[i] *= fb[i];
        fft(fa, true);

        res.resize(n);
        for (size_t i = 0; i < n; ++i)
            res[i] = int(fa[i].real() + 0.5);
    }
    normalize(res);
    for (int i = res.size() - 1; res[i] == 0 && i > 1; i--)
    {
        res.erase(res.begin() + i);
    }
}
void fft(std::vector<base>& a, bool invert) {
    int n = (int)a.size();
    if (n == 1)  return;

    std::vector<base> a0(n / 2), a1(n / 2);
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = a[i];
        a1[j] = a[i + 1];
    }
    fft(a0, invert);
    fft(a1, invert);

    double ang = 2 * PI / n * (invert ? -1 : 1);
    base w(1), wn(cos(ang), sin(ang));
    for (int i = 0; i < n / 2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (invert)
            a[i] /= 2, a[i + n / 2] /= 2;
        w *= wn;
    }
}
void normalize(std::vector<int>& l) {
    /*Нормализация числа - приведение каждого разряда в соответствие с системой счисления.
    *
    */
    for (int i = 0; i < l.size() - 1; ++i) {
        if (l[i] >= 10) { //если число больше максимального, то организовавается перенос
            int carryover = l[i] / 10;
            l[i + 1] += carryover;
            l[i] = l[i] % 10;
        }
        else if (l[i] < 0) { //если меньше - заем
            int carryover = (l[i] + 1) / 10 - 1;
            l[i + 1] += carryover;
            l[i] -= carryover * 10;
        }
    }
}