////////////////////////////////////////////////////////////////////////////////
/* 
 * File: drsa_bmp.c
 * 
 * \@brief Implementation of Dual Representation Simulated Annealing [1].
 * \@author Guilherme Oliveira Chagas (guilherme.o.chagas[a]gmail.com)
 * \@version 0.9
 * \@date This code was created on August 07, 2016, 02:03 PM
 * \@warning Apologizes about my bad english xD
 * \@copyright GNU General Public License
 * 
 * References:
 * [1] Torres-Jimenez et al. A dual representation simulated annealing algorithm
 * for the bandwidth minimization problem on graphs. Information Sciences, v.
 * 303, p. 33-49. 2015.
 * 
 * Disclaimer:
 * I am not a DRSA author, so it is possible that this DRSA version has errors 
 * and/or discrepancies with the actual Torres-Jimenez et al. DRSA algorithm.
 */
////////////////////////////////////////////////////////////////////////////////

// c library includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memcpy and memset
#include <time.h>
// file includes
#include "./mmio.h" // Matrix Market I/O library for ANSI c

/////////////////////////////////// Macros /////////////////////////////////////
// DRSA macros
#define PROBABILITY_REX 0.6 // probability of rex operator
#define PROBABILITY_NEX 0.2 // probability of nex operator
#define PROBABILITY_ROT 0.2 // probability of rot operator
#define T_F 1e-7 // final temperature of simulated annealing process
#define L 40 // initial lenght of Markov chain

// compressed sparse row (csr) matrix data structure
typedef struct csr
{
        double* nnz_; /**< non zero coefficients */
        int* col_ind_; /**< column indices */
        int* row_ind_; /**< row indices */
} csr;

// dual representation solution used in drsa algorithm (see [1])
typedef struct solution
{
        int* v_pi;
        int* v_rho;
} solution;

// global variables
int n = 0; // dimension of matrix

/////////////////////////////////// Methods ////////////////////////////////////

///////////////////////////////////////////////////
/* print csr row and col arrays indices */
///////////////////////////////////////////////////
void csr_print_indices(const csr* m_csr)
{
        int i, j;

        printf("Row indices: ");

        for (i = 0; i <= n; ++i)
                printf("%d, ", m_csr->row_ind_[i]);

        printf("\n");

        for (i = 0; i < n; ++i)
        {
                printf("Row %d: ", i);

                for(j = m_csr->row_ind_[i]; j < m_csr->row_ind_[i+1]; ++j)
                        printf("%d, ", m_csr->col_ind_[j]);

                printf("\n");
        }
}

///////////////////////////////////////////////////
/* print solution pi and rho arrays */
///////////////////////////////////////////////////
void print_solution(const solution* s)
{
        int i;
        printf("pi array: ");
        for (i = 0; i < n; ++i)
                printf("%d, ", s->v_pi[i]);
        printf("\nrho array: ");
        for (i = 0; i < n; ++i)
                printf("%d, ", s->v_rho[i]);
        printf("\n");
}

///////////////////////////////////////////////////
/* create a new copy solution from solution passed 
 * as parameter */
///////////////////////////////////////////////////
solution* copy_solution(const solution* x)
{
        solution* w = malloc(sizeof(*w));

        w->v_pi = (int*)calloc(n, sizeof(int));
        w->v_rho = (int*)calloc(n, sizeof(int));
        memcpy(w->v_pi, x->v_pi, n*sizeof(int));
        memcpy(w->v_rho, x->v_rho, n*sizeof(int));

        return w;
}

///////////////////////////////////////////////////
/* compute the matrix original bandwidth */
///////////////////////////////////////////////////
int compute_original_bandwidth(const csr* m_csr)
{
        int i, j;
        int b = 0;

        for (i = 0; i < n; ++i)
        {
                for (j = m_csr->row_ind_[i]; j < m_csr->row_ind_[i+1]; ++j)
                {
                        int b_i = abs(i - m_csr->col_ind_[j]);
                        if (b_i > b)
                                b = b_i;
                }
        }

        return b;
}

///////////////////////////////////////////////////
/* compute the matrix current bandwidth with 
 * solution v_pi */
///////////////////////////////////////////////////
int compute_current_bandwidth(const csr* m_csr, const int* v_pi)
{
        int i, j;
        int b = 0;

        for (i = 0; i < n; ++i)
        {
                for (j = m_csr->row_ind_[i]; j < m_csr->row_ind_[i+1]; ++j)
                {
                        int b_i = abs(v_pi[i] - v_pi[m_csr->col_ind_[j]]);
                        if(b_i > b)
                                b = b_i;		
                }
        }

        return b;
}

///////////////////////////////////////////////////
/* load matrix from file in matrix market format */
///////////////////////////////////////////////////
csr load_hb_matrix(const char* fname)
{
        MM_typecode matcode;
        FILE *f = NULL;
        csr csr_m;
        int i, j, m, nnz;
        int ret_code;

        if ((f = fopen(fname, "r")) == NULL)
        {
                exit(EXIT_FAILURE);
        }
 
        if (mm_read_banner(f, &matcode) != 0)
        {
                fprintf(stderr, "[ERROR] - could not process Matrix Market "
                        "banner");
                printf(" in file [%s]\n", fname);
                exit(EXIT_FAILURE);
        }

        // find out size of sparse matrix 
        if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nnz)) != 0)
        {
                fprintf(stderr, "ERROR - could not parse matrix size.\n");
                exit(EXIT_FAILURE);
        }

        if (m != n)
        {
                fprintf(stderr, "ERROR - matrix must be square! \n");
                exit(EXIT_FAILURE);
        }

        // initialize with zeros a mxn adjacencies matrix
        int** adj_matrix;
        adj_matrix = (int**)calloc(n, sizeof(int*));
        for (i = 0; i < n; ++i)
                adj_matrix[i] = (int *)calloc(m, sizeof(int));

        for (i = 0; i < nnz; ++i)
        {
                int	r = 0, // row index
                        c = 0; // column index
                double val = 0; // coefficient value
                if (mm_is_pattern(matcode)) // indices only
                {
                        if (fscanf(f, "%d %d \n", &r, &c) != 2)
                        {
                                fprintf(stderr, "ERROR - missing row and/or "
                                        "column data! \n");
                                exit(EXIT_FAILURE);
                        }
                }
                else
                {
                        if (fscanf(f, "%d %d %lg \n", &r, &c, &val) != 3)
                        {
                                fprintf(stderr, "ERROR - missing row and/or " 
                                        "column data! \n");
                                exit(EXIT_FAILURE);
                        }
                }
                --r; // adjust from 1-based to 0-based
                --c; //
                adj_matrix[r][c] = 1;
                if (mm_is_symmetric(matcode))
                        adj_matrix[c][r] = 1;
        }

        // parser to csr format
        csr_m.row_ind_ = (int *)calloc((n+1), sizeof(int));
        csr_m.col_ind_ = (int *)malloc((nnz) * sizeof(int));
        int n_r = 0;
        for (i = 0; i < n; ++i)
        {
                for (j = 0; j < n; ++j)
                {
                        if (adj_matrix[i][j] && i != j)
                                csr_m.col_ind_[n_r++] = j;
                }
                csr_m.row_ind_[i+1] = n_r;
        }

        // free memory
        free(f);
        for (i = 0; i < n; ++i)
        {
                free(adj_matrix[i]);
        }
        free(adj_matrix);

        return csr_m;
}

///////////////////////////////////////////////////
/* swap the idx1-th and idx2-th vertices labels */
///////////////////////////////////////////////////
void swap_labels(const solution* s, int idx1, int idx2)
{
        // swap rho labels
        int aux = s->v_rho[s->v_pi[idx1]];
        s->v_rho[s->v_pi[idx1]] = s->v_rho[s->v_pi[idx2]];
        s->v_rho[s->v_pi[idx2]] = aux;
        // swap pi labels
        aux = s->v_pi[idx1];
        s->v_pi[idx1] = s->v_pi[idx2];
        s->v_pi[idx2] = aux;
}

///////////////////////////////////////////////////
/* generate initial random solution for drsa */
///////////////////////////////////////////////////
solution* initial_random_solution()
{
        int i, j, k;
        solution* x = malloc(sizeof(solution));

        x->v_pi = (int*)calloc(n, sizeof(int));
        x->v_rho = (int*)calloc(n, sizeof(int));

        for (i = 1; i < n; ++i)
        {
                x->v_pi[i] = x->v_pi[i-1] + 1;
                x->v_rho[i] = x->v_rho[i-1] + 1;
        }

        // shuffle solution labels
        srand((double)time(NULL)); // seed for rand() function
        // at least n/2 swaps
        int number_of_changes = (n / 2 + (rand() % n)); 
        
        for (i = 0; i < number_of_changes; ++i)
        {
                j = (rand() % n);
                k = (rand() % n);

                if (j != k)
                {
                        swap_labels(x, j, k);
                }
        }

        for (i = 0; i < n; ++i)
        {
                if(i != x->v_rho[x->v_pi[i]])
                        printf("Insalubre!");
        }

        return x;
}

///////////////////////////////////////////////////
/* compute increment for Markov chain (see [1]) */
///////////////////////////////////////////////////
double compute_increment_factor(const double t_0, const double t_f, 
        const double l, const double l_f, const double alpha)
{
        double r = (log(t_f) - log(t_0)) / log(alpha);
        double gamma = exp((log(l_f) - log(l)) / r);

        return gamma;
}

///////////////////////////////////////////////////
/* random exchange operator (see [1]) */
///////////////////////////////////////////////////
void rex(const solution* w)
{
        int i = rand() % n;
        int j = rand() % n;
                        
        swap_labels(w, i, j);
}

///////////////////////////////////////////////////
/* neighborhood exchange operator (see [1]) */
///////////////////////////////////////////////////
void nex(const solution* w, const csr* m)
{
        unsigned i = rand() % n; // random vertex
        unsigned degree = m->row_ind_[i + 1] - m->row_ind_[i]; // number of adj
        unsigned j = rand() % degree; // random adjacent vertex

        // swap rho labels
        int aux = w->v_rho[w->v_pi[i]];
        w->v_rho[w->v_pi[i]] = w->v_rho[w->v_pi[m->col_ind_[m->row_ind_[i]+j]]];
        w->v_rho[w->v_pi[m->col_ind_[m->row_ind_[i] + j]]] = aux;

        // swap pi labels
        aux = w->v_pi[i];
        w->v_pi[i] = w->v_pi[m->col_ind_[m->row_ind_[i]+j]];
        w->v_pi[m->col_ind_[m->row_ind_[i]+j]] = aux;
}

///////////////////////////////////////////////////
/* rotation operator (see [1]) */
///////////////////////////////////////////////////
void rot(const solution* w)
{
        int k;
        int i = rand()%(n-5); // random vertex [0,n-5)
        int j = i + rand()%5 + 1; // random difference [1,5]

        int aux_rho = w->v_rho[i]; // change rho labels
        int aux_pi = w->v_pi[aux_rho]; // change pi labels
        for (k = i; k < j; ++k)
        {
                w->v_pi[w->v_rho[k]] = w->v_pi[w->v_rho[k+1]];
                w->v_rho[k] = w->v_rho[k+1];
        }

        w->v_pi[w->v_rho[j]] = aux_pi; 
        w->v_rho[j] = aux_rho;
}

///////////////////////////////////////////////////
/* neighborhood function g (see [1]) */
///////////////////////////////////////////////////
solution* neighborhood_function(const solution* x, const csr* m)
{
        double a = PROBABILITY_REX;
        double b = PROBABILITY_NEX;

        double p = (double)rand()/RAND_MAX; // random number [0,1]
        solution* w = copy_solution(x);
        
        if (p <= a)
                rex(w);
        else if (p <= a + b)
                nex(w,m);
        else
                rot(w);

        return w;
}

///////////////////////////////////////////////////
/* evaluation function f (see [1]) */
///////////////////////////////////////////////////
double evaluation_function(const csr* m, const solution* x, 
        const int*upsilon, int* d)
{
        int i, j;
        int b = 0;
        
        memset(d, 0, n*sizeof(int)); // reset array with 0

        for (i = 0; i < n; ++i)
        {
                for (j = m->row_ind_[i]; j < m->row_ind_[i+1]; ++j)
                {
                        int b_i = abs(x->v_pi[i] - x->v_pi[m->col_ind_[j]]);
                        ++d[b_i];
                        if(b_i > b)
                                b = b_i;
                }
        }

        double delta = 0;
        for (i = 0; i <= b; ++i)
                delta = (delta + d[i]) / upsilon[i];

        return (b + delta);
}

///////////////////////////////////////////////////
/* free memory pointers */
///////////////////////////////////////////////////
void free_solution(solution* s)
{
        free(s->v_pi);
        free(s->v_rho);
        free(s);
}

///////////////////////////////////////////////////
/* drsa algorithm (see [1]) */
///////////////////////////////////////////////////
void drsa(const csr* m, const double t_0, const double alpha, const double l_f)
{
        int i;
        
        double t_f = T_F;
        double l = L;
        double t = t_0;
        double gamma = compute_increment_factor(t_0, t_f, l, l_f, alpha);
        
        solution* x = initial_random_solution();
        solution* w = copy_solution(x);
        solution* y;

        // initialize arrays d and upsilon
        int* d = (int*)calloc(n, sizeof(int));
        int* upsilon = (int*)malloc(n*sizeof(int));
        for (i = 0; i < n; ++i)
                upsilon[i] = n+1-i; 

        double f_x = evaluation_function(m,x,upsilon,d);
        double f_w = f_x;
        double f_y;

        while (t > t_f)
        {
                int improvement = 0;
                for (i = 0; i < l; ++i)
                {
                        // flag to control free pointers
                        int x_point2_y = 0; 
                        y = neighborhood_function(x,m);
                        f_y = evaluation_function(m, y, upsilon,d);
                        
                        if (f_y < f_x)
                        {
                                // free pointers (for memory save)
                                free_solution(x); 
                                x = y;
                                x_point2_y = 1;
                                f_x = f_y;
                        }
                        
                        else if ((double)rand()/RAND_MAX < exp(-(f_y - f_x)/t))
                        {
                                // free pointers (for memory save)
                                free_solution(x); 
                                x = y;
                                x_point2_y = 1;
                                f_x = f_y;
                        }
                        
                        if (f_x < f_w)
                        {
                                // free pointers (for memory save)
                                free_solution(w); 
                                w = copy_solution(x); // save solution
                                f_w = f_x;
                                improvement = 1;
                        }
                        
                        if(!x_point2_y)
                                // free pointers (for memory save)
                                free_solution(y); 
                }

                if(!improvement)
                {
                        t *= alpha;
                        l *= gamma;
                }

                // printf("t = %f, l = %f, beta = %f\n", t, l, f_w);
        }
        printf("Best solution cost (beta) = %f\n", f_w);
}


///////////////////////////////////////////////////
/* checks if the input parameters are ok  */
///////////////////////////////////////////////////
int check_parameters(int argc, char **argv)
{
        if (argc != 3)
        {
                printf("[ERROR] Wrong number of parameters!\n");
                return 0;
        }

        if (strcmp(argv[1], "-f"))
        {
                printf("[ERROR] Wrong flag!\n");
                return 0;
        }

        FILE *f = NULL;
        if ((f = fopen(argv[2], "r")) == NULL) // checks if the path is valid
        {
                printf("[ERROR] Matrix file does not exists!\n");
                free(f);
                return 0;
        }

        free(f);
        return 1;
}

///////////////////////////////////////////////////
/* main routine */
///////////////////////////////////////////////////
int main(int argc, char **argv)
{
        if (!check_parameters(argc, argv))
        {
                printf("See README.md file...\n");
                return EXIT_FAILURE;
        }

        csr m = load_hb_matrix(argv[2]);

        printf ("Executing...\n");
        clock_t t = clock();

        // constants values following [1]
        const double alpha = 0.99;
        const double t_0 = 1000;
        const double l_f = 10 * n * m.row_ind_[n];

        drsa(&m, t_0, alpha, l_f);
        
        t = clock() - t;
        printf("It took me (%fs).\n", ((float)t)/CLOCKS_PER_SEC);

        return EXIT_SUCCESS;
}