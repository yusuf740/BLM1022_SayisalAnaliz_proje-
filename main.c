#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

static int is_space(char c)
{
    return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f';
}
static int is_digit(char c)
{
    return c >= '0' && c <= '9';
}
static int is_alpha(char c)
{
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

typedef enum
{
    T_NUMBER,
    T_VARIABLE,
    T_OP,
    T_FUNC,
    T_LPAREN,
    T_RPAREN
} TokenType;

typedef struct
{
    TokenType type;
    double value;
    char op;
    char *name;
} Token;

typedef struct Node
{
    TokenType type;
    double value;
    char op;
    char *name;
    struct Node *left, *right;
} Node;

Token *tokenize(const char *s, int *n_tokens);
void print_tokens(Token *toks, int n);
Node *parse_expression(Token *tokens, int *pos, int n);
double eval(Node *root, double x);
void free_tree(Node *n);
void free_tokens(Token *toks, int n);
Node *parse_primary(Token *toks, int *pos, int n);
Node *parse_factor(Token *toks, int *pos, int n);
Node *parse_term(Token *toks, int *pos, int n);
void bisection();
void regula_falsi();
void newton_raphson();
void find_inverse();
void cholesky();
void gauss_seidel_method();
double** take_matrix(int rows, int cols);
void free_matrix(double** matrix, int rows);
void clear_input_buffer();
double numerical_derivative(Node *func, double x, double h, int method);
void calculate_derivatives();
double simpson_1_3(Node *func, double a, double b, int n);
double simpson_3_8(Node *func, double a, double b, int n);
double trapezoidal_rule(Node *func, double a, double b, int n);
void calculate_integral();
void calculate_trapezoidal();
double factorial(int n);
double gregory_newton_interpolation(double* x_values, double* y_values, int n, double x);
void gregory_newton_method();

typedef struct {
    Node *tree;
    Token *tokens;
    int n_tokens;
} FunctionData;

void print_tokens(Token *toks, int n)
{
    for (int i = 0; i < n; i++)
    {
        Token *t = &toks[i];
        switch (t->type)
        {
        case T_NUMBER:
            printf("[NUM:%g] ", t->value);
            break;
        case T_VARIABLE:
            printf("[VAR:x] ");
            break;
        case T_OP:
            printf("[OP:%c] ", t->op);
            break;
        case T_FUNC:
            printf("[FUNC:%s] ", t->name);
            break;
        case T_LPAREN:
            printf("[ ( ] ");
            break;
        case T_RPAREN:
            printf("[ ) ] ");
            break;
        }
    }
    printf("\n");
}

FunctionData take_input_func(void)
{
    FunctionData result = {NULL, NULL, 0};
    size_t cap = 1024;
    char *buf = malloc(cap);
    if (!buf) {
        fprintf(stderr, "Memory allocation error\n");
        return result;
    }

    printf("Fonksiyonu girin (or: e^(x), sin(x*e^(5*x)), log_3(e^(x+2))):\n> ");
    fflush(stdout);

    if (!fgets(buf, cap, stdin)) {
        fprintf(stderr, "Error reading input\n");
        free(buf);
        return result;
    }
    size_t len = strlen(buf);
    if (len > 0 && buf[len-1] == '\n') buf[len-1] = '\0';
    if (buf[0] == '\0') {
        fprintf(stderr, "Empty input\n");
        free(buf);
        return result;
    }

    result.tokens = tokenize(buf, &result.n_tokens);
    free(buf);
    if (!result.tokens) return result;

    printf("=== Tokens ===\n");
    print_tokens(result.tokens, result.n_tokens);

    int pos = 0;
    result.tree = parse_expression(result.tokens, &pos, result.n_tokens);
    if (!result.tree || pos != result.n_tokens) {
        fprintf(stderr, "Parse error (pos=%d, tokens=%d)\n", pos, result.n_tokens);
        free_tree(result.tree);
        free_tokens(result.tokens, result.n_tokens);
        result.tree = NULL;
        result.tokens = NULL;
        result.n_tokens = 0;
    }
    return result;
}

int main(void)
{
    int choice = 0;
    int exit_flag = 0;
    size_t input_size = 256;
    char *input = (char *)malloc(input_size * sizeof(char));

    if (input == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    while (!exit_flag) {
        printf("\n1-Bisection\n");
        printf("2-Regula Falsi\n");
        printf("3-Newton Raphson\n");
        printf("4-Matrix Inverse\n");
        printf("5-Cholesky\n");
        printf("6-Gauss-Seidel\n");
        printf("7-Numerical Derivatives\n");
        printf("8-Simpson Integration\n");
        printf("9-Trapezoidal Integration\n");
        printf("10-Gregory-Newton Interpolation\n");
        printf("0-Exit\n");
        printf("Enter your choice (0-10): ");

        if (fgets(input, input_size, stdin) != NULL) {
            input[strcspn(input, "\n")] = 0;

            if (sscanf(input, "%d", &choice) != 1) {
                printf("Invalid input. Please enter a number.\n");
            } else {
                if (choice >=0 && choice <= 10){
                    switch (choice) {
                        case 1:
                            bisection();
                            clear_input_buffer();
                            break;
                        case 2:
                            regula_falsi();
                            clear_input_buffer();
                            break;
                        case 3:
                            newton_raphson();
                            clear_input_buffer();
                            break;
                        case 4:
                            find_inverse();
                            clear_input_buffer();
                            break;
                        case 5:
                            cholesky();
                            clear_input_buffer();
                            break;
                        case 6:
                            gauss_seidel_method();
                            clear_input_buffer();
                            break;
                        case 7:
                            calculate_derivatives();
                            clear_input_buffer();
                            break;
                        case 8:
                            calculate_integral();
                            clear_input_buffer();
                            break;
                        case 9:
                            calculate_trapezoidal();
                            clear_input_buffer();
                            break;
                        case 10:
                            gregory_newton_method();
                            clear_input_buffer();
                            break;
                        case 0:
                         printf("Exiting");
                            fflush(stdout);
                            for (int i = 0; i < 5; i++) {
                                for (volatile unsigned long j = 0; j < 100000000; j++);
                                printf(".");
                                fflush(stdout);
                            }
                            printf("\nGoodbye!\n");
                            exit_flag = 1;
                            break;
                    }
                } else {
                    printf("Invalid choice! Please enter a number between 0 and 10\n");
                }
            }
        } else {
            printf("Input error.\n");
            clearerr(stdin);
        }
    }

    free(input);
    return 0;
}

Token *tokenize(const char *s, int *n_tokens)
{
    const char *p = s;
    int cap = 256, nt = 0;
    Token *toks = malloc(cap * sizeof *toks);
    if (!toks) return NULL;

    while (*p) {
        if (is_space(*p)) {
            p++;
        }
        else if (*p == '-' && (nt == 0 || toks[nt-1].type == T_OP || toks[nt-1].type == T_LPAREN)) {
            if (nt + 1 > cap) { cap*=2; toks = realloc(toks, cap*sizeof *toks); }
            toks[nt++] = (Token){.type = T_OP, .op = 'u', .name = NULL};
            p++;
        }
        else if (*p == 'e' && *(p+1) == '^') {
            if (nt + 2 > cap) { cap*=2; toks = realloc(toks, cap*sizeof *toks); }
            toks[nt++] = (Token){.type = T_NUMBER, .value = M_E, .name = NULL};
            toks[nt++] = (Token){.type = T_OP, .op = '^', .name = NULL};
            p += 2;
        }
        else if (is_digit(*p) || (*p == '.' && is_digit(p[1]))) {
            char *end;
            double v = strtod(p, &end);
            if (nt >= cap) { cap*=2; toks = realloc(toks, cap*sizeof *toks); }
            toks[nt++] = (Token){.type = T_NUMBER, .value = v, .name = NULL};
            p = end;
        }
        else if (*p == 'x') {
            if (nt >= cap) { cap*=2; toks = realloc(toks, cap*sizeof *toks); }
            toks[nt++] = (Token){.type = T_VARIABLE, .name = NULL};
            p++;
        }
        else if (is_alpha(*p)) {
            const char *q = p;
            while (is_alpha(*q) || *q=='_' || is_digit(*q)) q++;
            int len = q - p;
            char *name = malloc(len+1);
            memcpy(name,p,len); name[len]='\0';
            if (nt >= cap) { cap*=2; toks = realloc(toks, cap*sizeof *toks); }
            toks[nt++] = (Token){.type = T_FUNC, .name = name, .value = 0};
            p = q;
        }
        else {
            Token tk = {0};
            if (strchr("+-*/^", *p)) { tk.type=T_OP; tk.op=*p; }
            else if (*p=='(') tk.type=T_LPAREN;
            else if (*p==')') tk.type=T_RPAREN;
            else {
                fprintf(stderr, "Unknown char '%c' at pos %ld\n", *p, p-s);
                free_tokens(toks, nt);
                return NULL;
            }
            if (nt >= cap) { cap*=2; toks = realloc(toks, cap*sizeof *toks); }
            toks[nt++] = tk;
            p++;
        }
    }
    *n_tokens = nt;
    return toks;
}

Node *new_node(Token tok)
{
    Node *n = malloc(sizeof *n);
    *n = (Node){.type=tok.type, .value=tok.value, .op=tok.op,
                .name= tok.name ? strdup(tok.name):NULL,
                .left=NULL, .right=NULL};
    return n;
}

Node *parse_expression(Token *toks,int *pos,int n)
{
    Node *left = parse_term(toks,pos,n);
    while (*pos<n && toks[*pos].type==T_OP && (toks[*pos].op=='+'||toks[*pos].op=='-')) {
        Token op = toks[(*pos)++];
        Node *right = parse_term(toks,pos,n);
        Node *p = new_node(op);
        p->left=left; p->right=right; left=p;
    }
    return left;
}

Node *parse_term(Token *toks,int *pos,int n)
{
    Node *left = parse_factor(toks,pos,n);
    while (*pos<n && toks[*pos].type==T_OP && (toks[*pos].op=='*'||toks[*pos].op=='/')) {
        Token op=toks[(*pos)++];
        Node *right=parse_factor(toks,pos,n);
        Node *p=new_node(op);
        p->left=left; p->right=right; left=p;
    }
    return left;
}

Node *parse_factor(Token *toks,int *pos,int n)
{
    Node *left = parse_primary(toks,pos,n);
    while (*pos<n && toks[*pos].type==T_OP && toks[*pos].op=='^') {
        Token op = toks[(*pos)++];
        Node *right=parse_factor(toks,pos,n);
        Node *p=new_node(op);
        p->left=left; p->right=right; left=p;
    }
    return left;
}

Node *parse_primary(Token *toks,int *pos,int n)
{
    if (*pos>=n) return NULL;
    Token tk=toks[*pos];
    if (tk.type==T_OP && tk.op=='u') {
        (*pos)++;
        Node *child = parse_primary(toks,pos,n);
        Node *p = new_node(tk);
        p->left = child;
        return p;
    }
    if (tk.type==T_NUMBER || tk.type==T_VARIABLE) {
        (*pos)++;
        return new_node(tk);
    }
    if (tk.type==T_FUNC) {
        (*pos)++;
        if (*pos<n && toks[*pos].type==T_LPAREN) {
            (*pos)++;
            Node *arg=parse_expression(toks,pos,n);
            if (*pos<n && toks[*pos].type==T_RPAREN) {
                (*pos)++;
                Node *fn=new_node(tk);
                fn->left=arg;
                return fn;
            }
        }
        return NULL;
    }
    if (tk.type==T_LPAREN) {
        (*pos)++;
        Node *e=parse_expression(toks,pos,n);
        if (*pos<n && toks[*pos].type==T_RPAREN) (*pos)++;
        return e;
    }
    return NULL;
}

double eval(Node *r, double x)
{
    double result = 0.0;
    
    if (!r) return result;
    
    if (r->type == T_NUMBER) {
        result = r->value;
    }
    else if (r->type == T_VARIABLE) {
        result = x;
    }
    else if (r->type == T_OP) {
        if (r->op == 'u') {
            result = -eval(r->left, x);
        }
        else {
            double L = eval(r->left, x), R = eval(r->right, x);
            if (r->op == '+') {
                result = L + R;
            }
            else if (r->op == '-') {
                result = L - R;
            }
            else if (r->op == '*') {
                result = L * R;
            }
            else if (r->op == '/') {
                if (fabs(R) < 1e-10) {
                    fprintf(stderr, "Error: Division by zero detected\n");
                    result = 0.0;
                } else {
                    result = L / R;
                }
            }
            else if (r->op == '^') {
                result = pow(L, R);
            }
        }
    }
    else if (r->type == T_FUNC) {
        double v = eval(r->left, x);
        if (!strcmp(r->name, "sin")) {
            result = sin(v);
        }
        else if (!strcmp(r->name, "cos")) {
            result = cos(v);
        }
        else if (!strcmp(r->name, "tan")) {
            result = tan(v);
        }
        else if (!strcmp(r->name, "arcsin") || !strcmp(r->name, "asin")) {
            if (v < -1.0 || v > 1.0) {
                fprintf(stderr, "Error: arcsin argument out of range [-1,1]: %f\n", v);
                result = 0.0;
            } else {
                result = asin(v);
            }
        }
        else if (!strcmp(r->name, "arccos") || !strcmp(r->name, "acos")) {
            if (v < -1.0 || v > 1.0) {
                fprintf(stderr, "Error: arccos argument out of range [-1,1]: %f\n", v);
                result = 0.0;
            } else {
                result = acos(v);
            }
        }
        else if (!strcmp(r->name, "arctan") || !strcmp(r->name, "atan")) {
            result = atan(v);
        }
        else if (!strcmp(r->name, "arccot") || !strcmp(r->name, "acot")) {
            if (fabs(v) < 1e-10) {
                fprintf(stderr, "Error: arccot argument too close to zero\n");
                result = 0.0;
            } else {
                result = atan(1.0 / v);
            }
        }
        else if (!strcmp(r->name, "exp")) {
            result = exp(v);
        }
        else if (!strcmp(r->name, "ln")) {
            if (v <= 0) {
                fprintf(stderr, "Error: ln argument must be positive (value: %.6f)\n", v);
                result = 0.0;
            } else {
                result = log(v);
            }
        }
        else if (!strcmp(r->name, "log")) {
            if (v <= 0) {
                fprintf(stderr, "Error: log argument must be positive (value: %.6f)\n", v);
                result = 0.0;
            } else {
                result = log10(v);
            }
        }
        else if (!strncmp(r->name, "log_", 4)) {
            if (r->name[4] == 'x' && r->name[5] == '\0') {
                double base = x;
                if (base <= 0 || base == 1) {
                    fprintf(stderr, "Error: Invalid logarithm base (x = %.6f)\n", base);
                    result = 0.0;
                }
                else {
                    result = log(v) / log(base);
                }
            } 
            else {
                double base = atof(r->name + 4);
                if (base <= 0 || base == 1) {
                    fprintf(stderr, "Error: Invalid logarithm base\n");
                    result = 0.0;
                }
                else {
                    result = log(v) / log(base);
                }
            }
        }
    }
    
    return result;
}

void free_tree(Node *n)
{
    if (!n) return;
    free_tree(n->left);
    free_tree(n->right);
    free(n->name);
    free(n);
}

void free_tokens(Token *toks,int n)
{
    for(int i=0;i<n;i++) free(toks[i].name);
    free(toks);
}

void clear_input_buffer() {
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF);
}

void bisection(){
    FunctionData func = take_input_func();
    if (!func.tree) {
        printf("Error: Invalid function input\n");
        return;
    }
    double x0, x1, expected_error;
    printf("Please enter x0 and x1: ");
    printf("\n(please enter it as x0 x1 with a space between them)\n"); 
    fflush(stdout);
    if (scanf("%lf %lf", &x0, &x1) != 2) {
        printf("Error: Invalid input for x0 and x1\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    printf("Please enter expected error: "); 
    fflush(stdout);
    if (scanf("%lf", &expected_error) != 1) {
        printf("Error: Invalid input for expected error\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    double f0 = eval(func.tree, x0);
    double f1 = eval(func.tree, x1);
    printf("f0: %lf, f1: %lf\n", f0, f1);
    if (f0 * f1 >= 0) { 
        fprintf(stderr, "Error: f(x0)*f(x1) >= 0\n"); 
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    double xm = 0.0, fm = 0.0;
    int iter = 0;
    const int MAX_ITER = 1000;
    int converged = 0; 
    
    printf("\nIteration\tx0\t\tx1\t\txm\t\tf(xm)\t\tError\n");
    printf("----------------------------------------------------------------\n");
    
    while (fabs(x1 - x0) > expected_error && iter < MAX_ITER && !converged) {
        xm = 0.5 * (x0 + x1);
        fm = eval(func.tree, xm);
        printf("%d\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               iter, x0, x1, xm, fm, fabs(x1 - x0));
               
        if (fabs(fm) < 1e-10) {
            converged = 1; 
        } else {
            if (f0 * fm < 0) {
                x1 = xm;
                f1 = fm;
            } else {
                x0 = xm;
                f0 = fm;
            }
        }
        iter++;
    }
    
    if (iter >= MAX_ITER) {
        fprintf(stderr, "Warning: Maximum iterations reached\n");
    }
    
    printf("Kok (bisection): %lf\n", xm);
    free_tree(func.tree);
    free_tokens(func.tokens, func.n_tokens);
}

void regula_falsi()
{
    FunctionData func = take_input_func();
    if (!func.tree) {
        printf("Error: Invalid function input\n");
        return;
    }

    double x0, x1, expected_error;
    printf("Please enter x0 and x1: "); 
    fflush(stdout);
    if (scanf("%lf %lf", &x0, &x1) != 2) {
        printf("Error: Invalid input for x0 and x1\n");
        printf("\n(please enter it as x0 x1 with a space between them)\n"); 
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }

    printf("Please enter expected error: "); 
    fflush(stdout);
    if (scanf("%lf", &expected_error) != 1) {
        printf("Error: Invalid input for expected error\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }

    double f0 = eval(func.tree, x0);
    double f1 = eval(func.tree, x1);
    printf("f0: %lf, f1: %lf\n", f0, f1);

    if (f0 * f1 >= 0) { 
        fprintf(stderr, "Error: f(x0)*f(x1) >= 0\n"); 
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }

    double xm = 0.0, fm = 0.0;
    int iter = 0;
    const int MAX_ITER = 1000;
    int converged = 0; 
    
    printf("\nIteration\tx0\t\tx1\t\txm\t\tf(xm)\t\tError\n");
    printf("----------------------------------------------------------------\n");
    
    while (fabs(x1 - x0) > expected_error && iter < MAX_ITER && !converged) {
        xm = (x1 * f0 - x0 * f1) / (f0 - f1);
        fm = eval(func.tree, xm);
        
        printf("%d\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               iter, x0, x1, xm, fm, fabs(x1 - x0));
        
        if (fabs(fm) < 1e-10) {
            converged = 1; 
        } else {
            if (f0 * fm < 0) {
                x1 = xm;
                f1 = fm;
            } else {
                x0 = xm;
                f0 = fm;
            }
        }
        iter++;
    }
    
    if (iter >= MAX_ITER) {
        fprintf(stderr, "Warning: Maximum iterations reached\n");
    }
    
    printf("Kok (regula falsi): %lf\n", xm);
    free_tree(func.tree);
    free_tokens(func.tokens, func.n_tokens);
}

void newton_raphson()
{
    FunctionData func = take_input_func();
    if (!func.tree) {
        printf("Error: Invalid function input\n");
        return;
    }

    double x0, expected_error;
    printf("Please enter x0: "); 
    fflush(stdout);
    if (scanf("%lf", &x0) != 1) {
        printf("Error: Invalid input for x0\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }

    printf("Please enter expected error: "); 
    fflush(stdout);
    if (scanf("%lf", &expected_error) != 1) {
        printf("Error: Invalid input for expected error\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    double f0 = eval(func.tree, x0);
    if (fabs(f0) < 1e-10) {
        printf("Root (Newton Raphson): %lf\n", x0);
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    double xm = 0.0, fm = 0.0;
    int iter = 0;
    const int MAX_ITER = 1000;
    const double h = 1e-6;
    int converged = 0; 
    
    printf("\nIteration\tx0\t\tf(x0)\t\tf'(x0)\t\txm\t\tf(xm)\t\tError\n");
    printf("----------------------------------------------------------------\n");
    
    while (fabs(f0) > expected_error && iter < MAX_ITER && !converged) {
        double f_plus = eval(func.tree, x0 + h);
        double f_minus = eval(func.tree, x0 - h);
        double df = (f_plus - f_minus) / (2 * h);
        
        if (fabs(df) < 1e-10) {
            fprintf(stderr, "Error: Derivative too close to zero\n");
            free_tree(func.tree);
            free_tokens(func.tokens, func.n_tokens);
            return;
        }
        
        xm = x0 - f0/df;
        fm = eval(func.tree, xm);
        
        printf("%d\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", iter, x0, f0, df, xm, fm, fabs(f0));
        
        if (fabs(fm) < 1e-10) {
            converged = 1;
        } else {
            x0 = xm;
            f0 = fm;
        }
        iter++;
    }
    
    if (iter >= MAX_ITER) {
        fprintf(stderr, "Warning: Maximum iterations reached\n");
    }
    
    printf("Kok (Newton Raphson): %lf\n", xm);
    free_tree(func.tree);
    free_tokens(func.tokens, func.n_tokens);
}

void find_inverse()
{    
    int rows, cols;
    printf("Enter the number of rows: ");
    if (scanf("%d", &rows) != 1 || rows <= 0) {
        printf("Error: Invalid number of rows\n");
        return;
    }

    printf("Enter the number of columns: ");
    if (scanf("%d", &cols) != 1 || cols <= 0) {
        printf("Error: Invalid number of columns\n");
        return;
    }

    if (rows != cols) {
        printf("Error: Matrix must be square for inverse calculation\n");
        return;
    }

    double** matrix = take_matrix(rows, cols);
    if (!matrix) return;

    double** aug = (double**)malloc(rows * sizeof(double*));
    if (!aug) {
        printf("Error: Memory allocation failed\n");
        free_matrix(matrix, rows);
        return;
    }

    for (int i = 0; i < rows; i++) {
        aug[i] = (double*)malloc(2 * cols * sizeof(double));
        if (!aug[i]) {
            printf("Error: Memory allocation failed\n");
            for (int j = 0; j < i; j++) {
                free(aug[j]);
            }
            free(aug);
            free_matrix(matrix, rows);
            return;
        }
        // Copy the original matrix
        for (int j = 0; j < cols; j++) {
            aug[i][j] = matrix[i][j];
        }
        // Add identity matrix
        for (int j = cols; j < 2 * cols; j++) {
            aug[i][j] = (j - cols == i) ? 1.0 : 0.0;
        }
    }

    printf("Augmented Matrix being processed:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 2 * cols; j++) {
            printf("%.6f ", aug[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < rows; i++) {
        int pivot = i;
        for (int j = i + 1; j < rows; j++) {
            if (fabs(aug[j][i]) > fabs(aug[pivot][i])) {
                pivot = j;
            }
        }

        if (fabs(aug[pivot][i]) < 1e-10) {
            printf("Matrix is singular (cannot be inverted)\n");
            free_matrix(aug, rows);
            free_matrix(matrix, rows);
            return;
        }

        if (pivot != i) {
            for (int j = 0; j < 2 * cols; j++) {
                double temp = aug[i][j];
                aug[i][j] = aug[pivot][j];
                aug[pivot][j] = temp;
            }
        }

        double pivot_val = aug[i][i];
        for (int j = 0; j < 2 * cols; j++) {
            aug[i][j] /= pivot_val;
        }

        for (int k = 0; k < rows; k++) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * cols; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }

    printf("\nTers Matris:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.6f ", aug[i][j + cols]);
        }
        printf("\n");
    }

    free_matrix(aug, rows);
    free_matrix(matrix, rows);
}

void cholesky()
{
    int n;
    printf("Enter matrix size: ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Error: Invalid matrix size\n");
        return;
    }
    
    double** A = take_matrix(n, n);
    if (!A) return;
    
    double* b = (double*)malloc(n * sizeof(double));
    if (!b) {
        printf("Memory allocation failed\n");
        free_matrix(A, n);
        return;
    }
    
    printf("Enter vector b:\n");
    for (int i = 0; i < n; i++) {
        printf("b[%d] = ", i);
        if (scanf("%lf", b + i) != 1) {
            printf("Error: Invalid input for b[%d]\n", i);
            free_matrix(A, n);
            free(b);
            return;
        }
    }
    
    double* x = (double*)calloc(n, sizeof(double));
    if (!x) {
        printf("Memory allocation failed\n");
        free_matrix(A, n);
        free(b);
        return;
    }
    
    double** L = (double**)malloc(n * sizeof(double*));
    double* y = (double*)malloc(n * sizeof(double));
    
    if (!L || !y) {
        printf("Memory allocation failed\n");
        if (L) free_matrix(L, n);
        if (y) free(y);
        free_matrix(A, n);
        free(b);
        free(x);
        return;
    }

    for (int i = 0; i < n; i++) {
        *(L + i) = (double*)calloc(n, sizeof(double));
        if (!*(L + i)) {
            printf("Memory allocation failed\n");
            for (int j = 0; j < i; j++) {
                free(*(L + j));
            }
            free(L);
            free(y);
            free_matrix(A, n);
            free(b);
            free(x);
            return;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += *(*(L + i) + k) * *(*(L + j) + k);
            }

            if (i == j) {
                double diag = *(*(A + i) + i) - sum;
                if (diag <= 0.0) {
                    printf("Error: Matrix is not positive definite\n");
                    free_matrix(L, n);
                    free(y);
                    free_matrix(A, n);
                    free(b);
                    free(x);
                    return;
                }
                *(*(L + i) + i) = sqrt(diag);
            } else {
                *(*(L + i) + j) = (*(*(A + i) + j) - sum) / *(*(L + j) + j);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += *(*(L + i) + k) * *(y + k);
        }
        *(y + i) = (*(b + i) - sum) / *(*(L + i) + i);
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int k = i + 1; k < n; k++) {
            sum += *(*(L + k) + i) * *(x + k);
        }
        *(x + i) = (*(y + i) - sum) / *(*(L + i) + i);
    }

    printf("\nSolution vector x:\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.6f\n", i, *(x + i));
    }

    free_matrix(L, n);
    free(y);
    free_matrix(A, n);
    free(b);
    free(x);
}

double** take_matrix(int rows, int cols)
{
    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) {
        printf("Memory allocation failed\n");
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        *(matrix + i) = (double*)malloc(cols * sizeof(double));
        if (!*(matrix + i)) {
            printf("Memory allocation failed\n");
            for (int j = 0; j < i; j++) {
                free(*(matrix + j));
            }
            free(matrix);
            return NULL;
        }
    }

    printf("Enter matrix elements:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("matrix[%d][%d] = ", i, j);
            if (scanf("%lf", *(matrix + i) + j) != 1) {
                printf("Error: Invalid input for matrix[%d][%d]\n", i, j);
                free_matrix(matrix, i);
                return NULL;
            }
        }
    }
    return matrix;
}

void free_matrix(double** matrix, int rows)
{
    if (matrix){
    for (int i = 0; i < rows; i++) {
        free(*(matrix + i));
    }
    free(matrix);
    }
    else {printf("Matrix is null\n");
}}

void gauss_seidel_method() {
    int n;
    printf("Enter matrix size: ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Error: Invalid matrix size\n");
        return;
    }
    
    double** A = take_matrix(n, n);
    if (!A) return;
    
    double* b = (double*)malloc(n * sizeof(double));
    if (!b) {
        printf("Memory allocation failed\n");
        free_matrix(A, n);
        return;
    }
    
    printf("Enter vector b:\n");
    for (int i = 0; i < n; i++) {
        printf("b[%d] = ", i);
        if (scanf("%lf", b + i) != 1) {
            printf("Error: Invalid input for b[%d]\n", i);
            free_matrix(A, n);
            free(b);
            return;
        }
    }
    
    double* x = (double*)calloc(n, sizeof(double));
    if (!x) {
        printf("Memory allocation failed\n");
        free_matrix(A, n);
        free(b);
        return;
    }
    
    int max_iter;
    printf("Enter maximum number of iterations: ");
    if (scanf("%d", &max_iter) != 1 || max_iter <= 0) {
        printf("Error: Invalid number of iterations\n");
        free_matrix(A, n);
        free(b);
        free(x);
        return;
    }
    
    double tol;
    printf("Enter convergence tolerance: ");
    if (scanf("%lf", &tol) != 1 || tol < 0) {
        printf("Error: Invalid tolerance\n");
        free_matrix(A, n);
        free(b);
        free(x);
        return;
    }
    
    int iter_count = 0;
    double* x_new = (double*)malloc(n * sizeof(double));
    if (!x_new) {
        printf("Memory allocation failed\n");
        free_matrix(A, n);
        free(b);
        free(x);
        return;
    }
    
    double err = tol + 1.0; 
    int i, j, iter = 0;
    

    for (i = 0; i < n; i++) {
        *(x_new + i) = *(x + i);
    }
    
    printf("\nIteration\tError\n");
    printf("-----------------------------------\n");

    while (iter < max_iter && err > tol) {

        for (i = 0; i < n; i++) {
            double sum1 = 0.0;  
            double sum2 = 0.0;  
            
            for (j = 0; j < i; j++) {
                sum1 += *(*(A + i) + j) * *(x_new + j);
            }
            
            for (j = i + 1; j < n; j++) {
                sum2 += *(*(A + i) + j) * *(x + j);
            }
            
            if (fabs(*(*(A + i) + i)) < 1e-10) {
                printf("Error: Diagonal element too close to zero. Method may not converge.\n");
                free(x_new);
                free_matrix(A, n);
                free(b);
                free(x);
                return;
            }
            
            *(x_new + i) = (*(b + i) - sum1 - sum2) / *(*(A + i) + i);
        }
        
        err = 0.0;
        for (i = 0; i < n; i++) {
            double diff = fabs(*(x_new + i) - *(x + i));
            if (diff > err) {
                err = diff;
            }
            *(x + i) = *(x_new + i);
        }
        
        iter++;
        iter_count = iter;
        printf("%d\t\t%.9f\n", iter, err);
    }
    
    if (iter >= max_iter) {
        printf("Warning: Maximum iterations reached without convergence\n");
    } else {
        printf("Solution converged in %d iterations\n", iter_count);
    }
    
    printf("\nSolution vector x:\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.9f\n", i, *(x + i));
    }
    
    free(x_new);
    free_matrix(A, n);
    free(b);
    free(x);
}

double numerical_derivative(Node *func, double x, double h, int method) {
    double result = 0.0;
    
    switch (method) {
        case 0: 
            result = (eval(func, x + h) - eval(func, x)) / h;
            break;
        case 1: 
            result = (eval(func, x) - eval(func, x - h)) / h;
            break;
        case 2: 
            result = (eval(func, x + h) - eval(func, x - h)) / (2 * h);
            break;
        default:
            printf("Invalid method for numerical differentiation\n");
            return 0;
    }
    
    return result;
}

void calculate_derivatives() {
    FunctionData func = take_input_func();
    if (!func.tree) {
        printf("Error: Invalid function input\n");
        return;
    }

    double x0, h;
    printf("Enter the point (x) at which to calculate derivative: ");
    fflush(stdout);
    if (scanf("%lf", &x0) != 1) {
        printf("Error: Invalid input for x\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }

    printf("Enter step size (h) for numerical differentiation: "); 
    fflush(stdout);
    if (scanf("%lf", &h) != 1 || h < 0) {
        printf("Error: Invalid input for h\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }

    double f0 = eval(func.tree, x0);
    
    printf("\nFunction value at x = %.6f: f(x) = %.9f\n\n", x0, f0);
    
    printf("First Derivatives:\n");
    printf("------------------\n");
    printf("Forward difference:  f'(x) = %.9f\n", 
           numerical_derivative(func.tree, x0, h, 0));
    printf("Backward difference: f'(x) = %.9f\n", 
           numerical_derivative(func.tree, x0, h, 1));
    printf("Central difference:  f'(x) = %.9f\n", 
           numerical_derivative(func.tree, x0, h, 2));
    
    free_tree(func.tree);
    free_tokens(func.tokens, func.n_tokens);
}

double simpson_1_3(Node *func, double a, double b, int n) {
    if (n % 2 != 0) {
        printf("Error: For Simpson 1/3 rule, n must be even.\n");
        return 0.0;
    }

    double h = (b - a) / n;
    double sum = eval(func, a) + eval(func, b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2 * eval(func, x);
        } else {
            sum += 4 * eval(func, x);
        }
    }
    return (h / 3.0) * sum;
}

double simpson_3_8(Node *func, double a, double b, int n) {
    if (n % 3 != 0) {
        printf("Error: For Simpson 3/8 rule, n must be divisible by 3.\n");
        return 0.0;
    }

    double h = (b - a) / n;
    double sum = eval(func, a) + eval(func, b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 3 == 0) {
            sum += 2 * eval(func, x);
        } else {
            sum += 3 * eval(func, x);
        }
    }

    return (3.0 * h / 8.0) * sum;
}

double trapezoidal_rule(Node *func, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = eval(func, a) + eval(func, b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += 2 * eval(func, x);
    }
    
    return (h / 2.0) * sum;
}

void calculate_integral() {
    FunctionData func = take_input_func();
    if (!func.tree) {
        printf("Error: Invalid function input\n");
        return;
    }

    double a, b;
    int n, method;
    
    printf("Enter the lower limit of integration (a): ");
    fflush(stdout);
    if (scanf("%lf", &a) != 1) {
        printf("Error: Invalid input for lower limit\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    printf("Enter the upper limit of integration (b): ");
    fflush(stdout);
    if (scanf("%lf", &b) != 1 || b <= a) {
        printf("Error: Invalid input for upper limit (must be greater than lower limit)\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    printf("Enter the number of intervals (n): ");
    fflush(stdout);
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Error: Invalid number of intervals\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    printf("Choose integration method:\n");
    printf("1 - Simpson's 1/3 rule (n must be even)\n");
    printf("2 - Simpson's 3/8 rule (n must be divisible by 3)\n");
    printf("Enter choice (1 or 2): ");
    fflush(stdout);
    if (scanf("%d", &method) != 1 || (method != 1 && method != 2)) {
        printf("Error: Invalid method choice\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    double result;
    if (method == 1) {
        if (n % 2 != 0) {
            printf("Error: For Simpson 1/3 rule, n must be even.\n");
            free_tree(func.tree);
            free_tokens(func.tokens, func.n_tokens);
            return;
        }
        result = simpson_1_3(func.tree, a, b, n);
        printf("\nSimpson's 1/3 Rule Result:\n");
    } 
    else {
        if (n % 3 != 0) {
            printf("Error: For Simpson 3/8 rule, n must be divisible by 3.\n");
            free_tree(func.tree);
            free_tokens(func.tokens, func.n_tokens);
            return;
        }
        result = simpson_3_8(func.tree, a, b, n);
        printf("\nSimpson's 3/8 Rule Result:\n");
    }
    
    printf("[%.6f to %.6f] f(x) dx = %.12f\n", a, b, result);
    
    free_tree(func.tree);
    free_tokens(func.tokens, func.n_tokens);
}

void calculate_trapezoidal() {
    FunctionData func = take_input_func();
    if (!func.tree) {
        printf("Error: Invalid function input\n");
        return;
    }

    double a, b;
    int n;
    
    printf("Enter the lower limit of integration (a): ");
    fflush(stdout);
    if (scanf("%lf", &a) != 1) {
        printf("Error: Invalid input for lower limit\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    printf("Enter the upper limit of integration (b): ");
    fflush(stdout);
    if (scanf("%lf", &b) != 1 || b <= a) {
        printf("Error: Invalid input for upper limit (must be greater than lower limit)\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    printf("Enter the number of intervals (n): ");
    fflush(stdout);
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Error: Invalid number of intervals\n");
        free_tree(func.tree);
        free_tokens(func.tokens, func.n_tokens);
        return;
    }
    
    double result = trapezoidal_rule(func.tree, a, b, n);
    printf("\nTrapezoidal Rule Result:\n");
    printf("[%.6f to %.6f] f(x) dx = %.12f\n", a, b, result);
    
    free_tree(func.tree);
    free_tokens(func.tokens, func.n_tokens);
}

double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

double gregory_newton_interpolation(double* x_values, double* y_values, int n, double x) {
    double h = x_values[1] - x_values[0];
    
    double** diff = (double**)malloc(n * sizeof(double*));
    if (!diff) {
        printf("Memory allocation failed\n");
        return 0.0;
    }
    
    for (int i = 0; i < n; i++) {
        diff[i] = (double*)malloc(n * sizeof(double));
        if (!diff[i]) {
            printf("Memory allocation failed\n");
            for (int j = 0; j < i; j++) {
                free(diff[j]);
            }
            free(diff);
            return 0.0;
        }
        diff[i][0] = y_values[i];
    }

    // Calculate forward differences table 
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
        }
    }


    double u = (x - x_values[0]) / h;
    double result = y_values[0];
    double u_term = 1.0;

    for (int i = 1; i < n; i++) {
        u_term *= (u - (i - 1));
        result += (u_term * diff[0][i]) / factorial(i);
    }

    for (int i = 0; i < n; i++) {
        free(diff[i]);
    }
    free(diff);

    return result;
}

void gregory_newton_method() {
    int n;
    printf("Enter the number of data points: ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Error: Invalid number of data points\n");
        return;
    }

    double* x_values = (double*)malloc(n * sizeof(double));
    double* y_values = (double*)malloc(n * sizeof(double));
    if (!x_values || !y_values) {
        printf("Memory allocation failed\n");
        free(x_values);
        free(y_values);
        return;
    }

    printf("Enter the data points (x y):\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        if (scanf("%lf", &x_values[i]) != 1) {
            printf("Error: Invalid input for x[%d]\n", i);
            free(x_values);
            free(y_values);
            return;
        }
        printf("y[%d] = ", i);
        if (scanf("%lf", &y_values[i]) != 1) {
            printf("Error: Invalid input for y[%d]\n", i);
            free(x_values);
            free(y_values);
            return;
        }
    }

    // Check if the x values are evenly spaced
    if (n > 2) {
        double h = x_values[1] - x_values[0];
        double tolerance = 1e-10;  
        int equal_spacing = 1;           
        for (int i = 2; i < n; i++) {
            double current_h = x_values[i] - x_values[i-1];
            if (fabs(current_h - h) > tolerance) {
                equal_spacing = 0;
                break;
            }
        }
        
        if (!equal_spacing) {
            printf("Warning: Gregory-Newton interpolation requires equally spaced x values.\n");
            printf("The data points you provided are not evenly spaced.\n");
            printf("Operation cancelled.\n");
            free(x_values);
            free(y_values);
            return;
        }
    }

    double x;
    printf("Enter the value of x for interpolation: ");
    if (scanf("%lf", &x) != 1) {
        printf("Error: Invalid input for x\n");
        free(x_values);
        free(y_values);
        return;
    }

    double result = gregory_newton_interpolation(x_values, y_values, n, x);
    printf("\nInterpolated value: f(%.6f) = %.9f\n", x, result);

    free(x_values);
    free(y_values);
}

