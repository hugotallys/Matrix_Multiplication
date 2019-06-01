#include <iostream>

/*
Strassen's Algorithm for matrix multiplication. To perform the matricial multiplication:

|1 2 3|     |4|   |35|
|4 5 6|  X  |5| = |83|
            |7|

Input must be like this:

2 3 -> number of lines and rows of the first matrix.
3 1 -> number of lines and rows of the second matrix.

1 2 3 4 5 6 -> elements of the first matrix (line by line).
4 5 7 -> elements of the second matrix (line by line).
*/

// Allocates a nxn square matrix filled with zeroes.
int * allocate_matrix(unsigned n)
{
    unsigned i;
    int *zeroes = new int[n*n];

    for(i = 0; i < n*n; ++i)
        zeroes[i] = 0;
    
    return zeroes;
}

// Takes a pointer to a square matrix and deallocates it.
void delete_matrix(int *matrix)
{
    delete[] matrix;
}

// Takes a square matrix nxn and returns the minor matrix lxc.
int * trim_matrix(int *matrix, unsigned n, unsigned l, unsigned c)
{
    int *t;
    unsigned i , j;

    t = new int[l*c];

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
void print_matrix(int *matrix, unsigned l, unsigned c)
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
int * add_matrix(int *a, int *b, unsigned n)
{
    int *s;
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
int * sub_matrix(int *a, int *b, unsigned n)
{
    int *s;
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
// Standard algorithm (O(n^3)) for matrix multiplication with a divide-and-conquer approach.
int * mult_matrix(int *a, int *b, unsigned n)
{
    int *c = allocate_matrix(n);

    if(n == 1)
    {
       c[0] = a[0] * b[0];
    }
    else
    {
        unsigned i, j, m;

        m = n >> 1;

        int *a11 = allocate_matrix(m);
        int *a12 = allocate_matrix(m);
        int *a21 = allocate_matrix(m);
        int *a22 = allocate_matrix(m);

        int *b11 = allocate_matrix(m);
        int *b12 = allocate_matrix(m);
        int *b21 = allocate_matrix(m);
        int *b22 = allocate_matrix(m);

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

        int *c11 = add_matrix(mult_matrix(a11, b11, m), mult_matrix(a12, b21, m), m);
        int *c12 = add_matrix(mult_matrix(a11, b12, m), mult_matrix(a12, b22, m), m);
        int *c21 = add_matrix(mult_matrix(a21, b11, m), mult_matrix(a22, b21, m), m);
        int *c22 = add_matrix(mult_matrix(a21, b12, m), mult_matrix(a22, b22, m), m);

        delete_matrix(a11);
        delete_matrix(a12);
        delete_matrix(a21);
        delete_matrix(a22);
        delete_matrix(b11);
        delete_matrix(b12);
        delete_matrix(b22);
        delete_matrix(b22);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                c[i*n + j] = c11[i*m + j];
                c[i*n + (j+m)] = c12[i*m + j];
                c[(i+m)*n + j] = c21[i*m + j];
                c[(i+m)*n + (j+m)] = c22[i*m + j];
            }
        }
        
        delete_matrix(c11);
        delete_matrix(c12);
        delete_matrix(c21);
        delete_matrix(c22);
    }
    return c;
}

// Takes two square matrices and their size and returns their multiplication.
// Strassen's algorithm (O(n^(lg7))) for matrix multiplication.
int * strassen_algorithm(int *a, int *b, unsigned n)
{
    int *s = allocate_matrix(n);

    if(n == 1)
    {
        s[0] = a[0] * b[0];
    }
    else
    {
        unsigned i, j, m;

        m = n >> 1;

        int *a11 = allocate_matrix(m);
        int *a12 = allocate_matrix(m);
        int *a21 = allocate_matrix(m);
        int *a22 = allocate_matrix(m);

        int *b11 = allocate_matrix(m);
        int *b12 = allocate_matrix(m);
        int *b21 = allocate_matrix(m);
        int *b22 = allocate_matrix(m);

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

        int *p1 = strassen_algorithm(a11, sub_matrix(b12, b22, m), m);
        int *p2 = strassen_algorithm(add_matrix(a11, a12, m), b22, m);
        int *p3 = strassen_algorithm(add_matrix(a21, a22, m), b11, m);
        int *p4 = strassen_algorithm(a22, sub_matrix(b21, b11, m), m);
        int *p5 = strassen_algorithm(add_matrix(a11, a22, m), add_matrix(b11, b22, m), m);
        int *p6 = strassen_algorithm(sub_matrix(a12, a22, m), add_matrix(b21, b22, m), m);
        int *p7 = strassen_algorithm(sub_matrix(a11, a21, m), add_matrix(b11, b12, m), m);

        delete_matrix(a11);
        delete_matrix(a12);
        delete_matrix(a21);
        delete_matrix(a22);
        delete_matrix(b11);
        delete_matrix(b12);
        delete_matrix(b22);
        delete_matrix(b22);

        int *s11 = add_matrix(add_matrix(sub_matrix(p4, p2, m), p5, m), p6, m);
        int *s12 = add_matrix(p1, p2, m);
        int *s21 = add_matrix(p3, p4, m);
        int *s22 = add_matrix(sub_matrix(sub_matrix(p5, p3, m), p7, m), p1, m);

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

// Takes to positive numbers returns the value of the greater.
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
    int *A, *B, *S, *C;
    unsigned lin_a, col_a, lin_b, col_b, i, j, n, m;

    scanf("%u %u", &lin_a, &col_a);
    scanf("%u %u", &lin_b, &col_b);

    if(col_a != lin_b)
    {
        printf("Impossible to perform multiplication. Recall that the number of columns of the first matrix must be equals to the number of lines of the second matrix.\n");
        return 0;
    }

    n = max(lin_b, max(col_b, max(lin_a, col_a)));
    m = closest_power_of_two(n);

    A = allocate_matrix(m);
    B = allocate_matrix(m);

    for (i = 0; i < lin_a; i++)
    {
        for (j = 0; j < col_a; j++)
        {
            scanf("%d", A + (i*m) + j);
        }
    }

    for (i = 0; i < lin_b; i++)
    {
        for (j = 0; j < col_b; j++)
        {
            scanf("%d", B + (i*m) + j);
        }
    }

    S = strassen_algorithm(A, B, m);
    C = trim_matrix(S, m, lin_a, col_b);
    print_matrix(C, lin_a, col_b);

    return 0;
}