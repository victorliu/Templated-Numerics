INDEX i, j, k;
INDEX n1, n2;
INDEX ldf, ldg;
int TransF, TransG;
const BASE *F, *G;

if(alpha == 0.0 && beta == 1.0)
return;

if(Order == CblasRowMajor){
n1 = M;
n2 = N;
F = A;
ldf = lda;
TransF = (TransA == CblasConjTrans) ? CblasTrans : TransA;
G = B;
ldg = ldb;
TransG = (TransB == CblasConjTrans) ? CblasTrans : TransB;
}else{
n1 = N;
n2 = M;
F = B;
ldf = ldb;
TransF = (TransB == CblasConjTrans) ? CblasTrans : TransB;
G = A;
ldg = lda;
TransG = (TransA == CblasConjTrans) ? CblasTrans : TransA;
}

/* form  y := beta*y */
if(beta == 0.0){
for (i = 0; i < n1; i++){
for (j = 0; j < n2; j++){
C[ldc * i + j] = 0.0;
}
}
}else if(beta != 1.0){
for (i = 0; i < n1; i++){
for (j = 0; j < n2; j++){
C[ldc * i + j] *= beta;
}
}
}

if(alpha == 0.0)
return;

if(TransF == CblasNoTrans && TransG == CblasNoTrans){

/* form  C := alpha*A*B + C */

for (k = 0; k < K; k++){
for (i = 0; i < n1; i++){
const BASE temp = alpha * F[ldf * i + k];
if(temp != 0.0){
for (j = 0; j < n2; j++){
C[ldc * i + j] += temp * G[ldg * k + j];
}
}
}
}

}else if(TransF == CblasNoTrans && TransG == CblasTrans){

/* form  C := alpha*A*B' + C */

for (i = 0; i < n1; i++){
for (j = 0; j < n2; j++){
BASE temp = 0.0;
for (k = 0; k < K; k++){
temp += F[ldf * i + k] * G[ldg * j + k];
}
C[ldc * i + j] += alpha * temp;
}
}

}else if(TransF == CblasTrans && TransG == CblasNoTrans){

for (k = 0; k < K; k++){
for (i = 0; i < n1; i++){
const BASE temp = alpha * F[ldf * k + i];
if(temp != 0.0){
for (j = 0; j < n2; j++){
C[ldc * i + j] += temp * G[ldg * k + j];
}
}
}
}

}else if(TransF == CblasTrans && TransG == CblasTrans){

for (i = 0; i < n1; i++){
for (j = 0; j < n2; j++){
BASE temp = 0.0;
for (k = 0; k < K; k++){
temp += F[ldf * k + i] * G[ldg * j + k];
}
C[ldc * i + j] += alpha * temp;
}
}

}else {
BLAS_ERROR("unrecognized operation");
}
