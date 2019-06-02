#include <iostream>

// Allocates a nxn square matrix filled with zeroes.
long long int * allocate_matrix(unsigned n)
{
    unsigned i;
    long long int *zeroes = new long long int[n*n];

    for(i = 0; i < n*n; ++i)
        zeroes[i] = 0;
    
    return zeroes;
}

// Takes a pointer to a square matrix and deallocates it.
void delete_matrix(long long int *matrix)
{
    delete[] matrix;
}

// Takes a square matrix nxn and returns the minor matrix lxc.
long long int * trim_matrix(long long int *matrix, unsigned n, unsigned l, unsigned c)
{
    long long int *t;
    unsigned i , j;

    t = new long long int[l*c];

    for(i = 0; i < l; ++i)
    {
        for(j = 0; j < c; j++)
        {
            t[i*c + j] = matrix[i*n + j];
        }
    }

    delete_matrix(matrix);

    return t;
}

// Takes a square matrix and its dimensions and outputs the its content to the screen.
void print_matrix(long long int *matrix, unsigned l, unsigned c)
{
    unsigned i, j;
    
    for(i = 0; i < l; ++i)
    {
        for(j = 0; j < c; j++)
        {
            printf("%d ", matrix[i*c + j]);
        }
        printf("\n");
    }
}

// Takes two square matrices and their size and returns their sum.
long long int * add_matrix(long long int *a, long long int *b, unsigned n)
{
    long long int *s;
    unsigned i, j;

    s = allocate_matrix(n);
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < n; j++)
        {
            s[i*n + j] = a[i*n + j] + b[i*n + j];
        }
    }
    return s;
}

// Takes two square matrices and their size and returns their subtraction.
long long int * sub_matrix(long long int *a, long long int *b, unsigned n)
{
    long long int *s;
    unsigned i, j;

    s = allocate_matrix(n);
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < n; j++)
        {
            s[i*n + j] = a[i*n + j] - b[i*n + j];
        }
    }
    return s;
}

// Takes two square matrices and their size and returns their multiplication.
// Strassen's algorithm (O(n^(lg7))) for matrix multiplication.
long long int * strassen_algorithm(long long int *a, long long int *b, unsigned n)
{
    long long int *s = allocate_matrix(n);

    if(n == 1)
    {
        s[0] = a[0] * b[0];
    }
    else
    {
        unsigned i, j, m;

        m = n >> 1;

        long long int *a11 = allocate_matrix(m);
        long long int *a12 = allocate_matrix(m);
        long long int *a21 = allocate_matrix(m);
        long long int *a22 = allocate_matrix(m);

        long long int *b11 = allocate_matrix(m);
        long long int *b12 = allocate_matrix(m);
        long long int *b21 = allocate_matrix(m);
        long long int *b22 = allocate_matrix(m);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                a11[i*m + j] = a[i*n + j];
                a12[i*m + j] = a[i*n + j+m];
                a21[i*m + j] = a[(i+m)*n + j];
                a22[i*m + j] = a[(i + m)*n + (j+m)];
                b11[i*m + j] = b[i*n + j];
                b12[i*m + j] = b[i*n + j+m];
                b21[i*m + j] = b[(i+m)*n + j];
                b22[i*m + j] = b[(i + m)*n + (j+m)];
            }
        }

        long long int *p1 = strassen_algorithm(a11, sub_matrix(b12, b22, m), m);
        long long int *p2 = strassen_algorithm(add_matrix(a11, a12, m), b22, m);
        long long int *p3 = strassen_algorithm(add_matrix(a21, a22, m), b11, m);
        long long int *p4 = strassen_algorithm(a22, sub_matrix(b21, b11, m), m);
        long long int *p5 = strassen_algorithm(add_matrix(a11, a22, m), add_matrix(b11, b22, m), m);
        long long int *p6 = strassen_algorithm(sub_matrix(a12, a22, m), add_matrix(b21, b22, m), m);
        long long int *p7 = strassen_algorithm(sub_matrix(a11, a21, m), add_matrix(b11, b12, m), m);

        delete_matrix(a11);
        delete_matrix(a12);
        delete_matrix(a21);
        delete_matrix(a22);
        delete_matrix(b11);
        delete_matrix(b12);
        delete_matrix(b21);
        delete_matrix(b22);

        long long int *s11 = add_matrix(add_matrix(sub_matrix(p4, p2, m), p5, m), p6, m);
        long long int *s12 = add_matrix(p1, p2, m);
        long long int *s21 = add_matrix(p3, p4, m);
        long long int *s22 = add_matrix(sub_matrix(sub_matrix(p5, p3, m), p7, m), p1, m);

        delete_matrix(p1);
        delete_matrix(p2);
        delete_matrix(p3);
        delete_matrix(p4);
        delete_matrix(p5);
        delete_matrix(p6);
        delete_matrix(p7);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                s[i*n + j] = s11[i*m + j];
                s[i*n + (j+m)] = s12[i*m + j];
                s[(i+m)*n + j] = s21[i*m + j];
                s[(i+m)*n + (j+m)] = s22[i*m + j];
            }
        }

        delete_matrix(s11);
        delete_matrix(s12);
        delete_matrix(s21);
        delete_matrix(s22);
    }
    return s;    
}

// Takes two positive numbers returns the value of the greater.
unsigned max(unsigned a, unsigned b)
{
    return a > b ? a : b;
}

// Takes a postitve number n and returns the first power of two b such that b >= n. 
unsigned closest_power_of_two(unsigned n)
{
    unsigned b = 2;
    while(n > b)
        b = b << 1;
    return b;
}

int main()
{
    long long int *A, *B, *S, *C;
    unsigned lin_a, col_a, lin_b, col_b, i, j, n, m;

    scanf("%u %u %u %u", &lin_a, &col_a, &lin_b, &col_b);

    n = max(lin_b, max(col_b, max(lin_a, col_a)));
    m = closest_power_of_two(n);

    A = allocate_matrix(m);
    B = allocate_matrix(m);

    for (i = 0; i < lin_a; i++)
    {
        for (j = 0; j < col_a; j++)
        {
            scanf("%lld", A + (i*m) + j);
        }
    }

    for (i = 0; i < lin_b; i++)
    {
        for (j = 0; j < col_b; j++)
        {
            scanf("%lld", B + (i*m) + j);
        }
    }

    S = strassen_algorithm(A, B, m);
    C = trim_matrix(S, m, lin_a, col_b);
    printf("17111663\n%u %u\n", lin_a, col_b);
    print_matrix(C, lin_a, col_b);
    
    return 0;
}
