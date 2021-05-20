
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "omp.h"

using namespace std;

int main(int argc, char* argv[]);

int congruence(int a, int b, int c, bool* error);
int GCD(int i, int j);
int MAX(int i, int j);
int MIN(int i, int j);
int SIGN(int i);
int LCRG_EVALUATE(int a, int b, int c, int x);
void LCRG_ANBN(int a, int b, int c, int n, int* an, int* bn);
int POWER_MOD(int a, int n, int m);

int main(int argc, char* argv[]) {
    // Дана програма показує, як P-кількість процесорів може згенерувати
    // таку саму послідовність як і один процесор
    int a, an, b, bn, c, u, v;
    double Start, Finish;
    int const p = omp_get_num_procs();
    int *xi = new int[p];
    double duration = 0;

    cout << "    U = ( A * V + B ) mod C\n";//Правило, за яким шукається послідовність псевдовипадкових чисел
    cout << "\n";

    //параметри правила-функції
    a = 16807;
    b = 0;
    c = 2147483637;

        cout << "  LCRG parameters:\n";
        cout << "\n";
        cout << "  A  = " << a << "\n";
        cout << "  B  = " << b << "\n";
        cout << "  C  = " << c << "\n";

    //кількість згенерованих чисел
#define k_hi 40000000
    int *result = new int[k_hi];
    // Р-процесор обчислює параметри для Р-ої частини послідовності
    //початкове значення для запуску послідовності
    v = 12345;
    result[0] = v;
    LCRG_ANBN(a, b, c, p, &an, &bn);

    u = v;
    int q = 0;
    int w = 0;
    Start = omp_get_wtime();

#pragma omp parallel for  shared(xi) private(q,v)
    for (int k = 1; k < k_hi;k+=(k_hi/p)) {
        LCRG_ANBN(a, b, c, k, &an, &bn);                 
        v = LCRG_EVALUATE(an, bn, c, u);
        q = omp_get_num_threads();
        xi[q] = v;
        //cout << v << "\n";
    }
    
#pragma omp parallel for shared(xi,result) private(q,w,v)
    for (int k = 1; k < k_hi;k += (k_hi / p)) {
        q = omp_get_num_threads();
        u = xi[q];
        for (int i = 1; i < k_hi;i +=p) {                   
            v = LCRG_EVALUATE(a, b, c, u);
            u = v;
            //result[(k_hi * q + i )/p] = v;
            //cout << v <<"\n";
        } 
    }       
    
    Finish = omp_get_wtime();           

    duration = Finish - Start;

    cout << "\n";
    cout << "Time of parallel execution: " << duration;//час виконання паралельного коду

    return 0;
}

//Розв'язання рівняння ( A * X ) mod B = C
int congruence(int a, int b, int c, bool* error) {
#define N_MAX 100
    int a_copy;
    int a_mag;
    int a_sign;
    int b_copy;
    int b_mag;
    int c_copy;
    int g;
    int k;
    int n;
    int q[N_MAX];
    bool swap;
    int x, y, z;

    *error = false;
    x = 0;
    y = 0;
    
    if (a == 0 && b == 0 && c == 0) {
        x = 0;
        return x;
    }
    else if (a == 0 && b == 0 && c != 0) {
        *error = true;
        x = 0;
        return x;
    }
    else if (a == 0 && b != 0 && c == 0) {
        x = 0;
        return x;
    }
    else if (a == 0 && b != 0 && c != 0) {
        x = 0;
        if ((c % b) != 0) {
            *error = true;
        }
        return x;
    }
    else if (a != 0 && b == 0 && c == 0) {
        x = 0;
        return 0;
    }
    else if (a != 0 && b == 0 && c != 0) {
        x = c / a;
        if ((c % a) != 0) {
            *error = true;
        }
        return x;
    }
    else if (a != 0 && b != 0 && c == 0) {
        x = 0;
        return x;
    }
    //1. Обчислення НСД чисел a,b, які мають ділитися на С
    g = GCD(a, b);

    if ((c % g) != 0) {
        *error = true;
        return x;
    }
    a_copy = a / g;
    b_copy = b / g;
    c_copy = c / g;
    //2. Розділимо числа A,B на значення за модулем і знак
    a_mag = abs(a_copy);
    a_sign = SIGN(a_copy);
    b_mag = abs(b_copy);
    //перевірка чи a==1 або b==1
    if (a_mag == 1) {
        x = a_sign * c_copy;
        return x;
    }
    else if (b_mag == 1) {
        x = 0;
        return x;
    }
    //Обчислення послідовності Евклідової остачі
    if (b_mag <= a_mag) {
        swap = false;
        q[0] = a_mag;
        q[1] = b_mag;
    }
    else {
        swap = true;
        q[0] = b_mag;
        q[1] = a_mag;
    }
    n = 3;

    for (;;) {
        q[n - 1] = (q[n - 3] % q[n - 2]);

        if (q[n - 1] == 1) {
            break;
        }

        n += 1;
        if (N_MAX < n) {
            *error = true;
            cout << "\n";
            cout << "CONGRUENCE - Fatal error!\n";
            cout << "  Exceeded number of iterations.\n";
            exit(1);
        }
    }

    //4. Розв'язання рівняння X*A_MAG + Y*B_MAG = 1
    y = 0;
    for (k = n; k >= 2;k--) {
        x = y;
        y = (1 - x * q[k - 2]) / q[k - 1];
    }

    //5. Перевіряємо чи треба обмінювати значення між собою
    if (swap) {
        z = x;
        x = y;
        y = z;
    }
    //6. Повертаємо знак, щоб ми мали рівняння X * A + Y * B = 1
    x = x * a_sign;
    //7. Множимо на С, щоб отримати X * A + Y * B = C
    x = x * c_copy;
    //8. Беремо значення за модулем В
    x = x % b;
    //9. Якщо корінь від'ємний, робимо його дадатнім
    if (x < 0) {
        x += b;
    }

    return x;
#undef N_MAX;
}

//Знаходить НСД чисел i , j
int GCD(int i, int j) {
    int ip;
    int iq;
    int ir;
    //перевіряємо чи числа не є рівними 0
    if (i == 0) {
        return MAX(1, abs(j));
    }
    else if (j == 0) {
        return MAX(1, abs(i));
    }
    //Обираємо менше і більше числа
    ip = MAX(abs(i), abs(j));
    iq = MIN(abs(i), abs(j));
    //Запускаємо алгоритм Евкліда
    for (;;) {
        ir = ip % iq;
        if (ir == 0) {
            break;
        }
        ip = iq;
        iq = ir;
    }

    return iq;
}
//Шукає максимальне значення з двох
int MAX(int i, int j) {
    int value;
    if (j < i) {
        value = i;
    }
    else {
        value = j;
    }
    return value;
}
//Шукає мінімальне значення з двох
int MIN(int i, int j) {
    int value;
    if (i < j) {
        value = i;
    }
    else {
        value = j;
    }
    return value;
}
//Шукає знак обраного числа
int SIGN(int i) {
    int value;
    if (i < 0) {
        value = -1;
    }
    else if (i > 0) {
        value = 1;
    }
    else {
        value = 0;
    }
    return value;
}
//Обчислення ( A^N ) mod M 
//за допомогою бінарного піднесення до степеня по модулю
int POWER_MOD(int a, int n, int m) {
    long long int aa; //square
    int d;
    long long int m2;
    int x;
    long long int x2;

    if (a < 0) {
        return -1;
    }

    if (m <= 0) {
        return -1;
    }
    if (n < 0) {
        return -1;
    }
    //Означкємо квадрат числа самим числом 
    aa = (long long int) a;
    m2 = (long long int) m;
    x2 = (long long int) 1;

    while (n > 0) {
        d = n % 2;
        if (d == 1) {
            x2 = (x2 * aa) % m2;
        }

        aa = (aa * aa) % m2;
        n = (n - d) / 2;
    }
    //Переконаємося, що х більше рівне 0
    if (x2 < 0) {
        x2 += m2;
    }

    x = (int)x2;

    return x;
}
//Обчислює "N-степінь" функції генератора псевдовипадкових чисел
//F(X_i) = (a * F(X_i-1) + b) mod c
//F(N) = a^N * F(0) + ( a^N - 1) / ( a - 1 ) * b
//     = AN * F(0) + BN
void LCRG_ANBN(int a, int b, int c, int n, int* an, int* bn) {
    int am1;
    int anm1tb;
    bool ierror;

    if (n < 0) {
        cerr << "\n";
        cerr << "LCRG_ANBN - Fatal error!\n";
        cerr << "Illegal input value of N = " << n << "\n";
        exit(1);
    }

    if (c <= 0) {
        cerr << "\n";
        cerr << "LCRG_ANBN - Fatal error!\n";
        cerr << "Illegal input value of C = " << c << "\n";
        exit(1);
    }

    if (n == 0) {
        *an = 1;
        *bn = 0;
    }
    else if (n == 1) {
        *an = a;
        *bn = b;
    }
    else
    {
        //Обчислимо A^N
        *an = POWER_MOD(a, n, c);
        //Розв'яжемо (a - 1)* BN = (a ^ N - 1) mod B відносно for BN
        am1 = a - 1;
        anm1tb = (*an - 1)*b;
        //Використаємо метод для розв'язання модульних рівнянь
        *bn = congruence(am1, c, anm1tb, &ierror);

    }
    return;
}
//Обчислює саму функцію генератора y = ( A * x + B ) mod C за усіма значеннями 
int LCRG_EVALUATE(int a, int b, int c, int x) {
    long long int aLong;
    long long int bLong;
    long long int cLong;
    long long int xLong;
    int y;
    long long int yLong;

    aLong = (long long int) a;
    bLong = (long long int) b;
    cLong = (long long int) c;
    xLong = (long long int) x;

    yLong = (aLong * xLong + bLong) % cLong;
    y = (int)yLong;

    if (y < 0) {
        y += c;
    }

    return y;
}
