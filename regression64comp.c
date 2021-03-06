/**************************************************************************************
*
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
*
* Progetto dell'algoritmo di Regressione
* in linguaggio assembly x86-32 + SSE
*
* Fabrizio Angiulli, aprile 2019
*
**************************************************************************************/

/*
*
* Software necessario per l'esecuzione:
*
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
*
* entrambi sono disponibili come pacchetti software
* installabili mediante il packaging tool del sistema
* operativo; per esempio, su Ubuntu, mediante i comandi:
*
*    sudo apt-get install nasm
*    sudo apt-get install gcc
*
* potrebbe essere necessario installare le seguenti librerie:
*
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
*
* Per generare il file eseguibile:
*
* nasm -f elf32 regression32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o regression32.o regression32c.c -o regression32c -lm && ./regression32c $pars
*
* oppure
*
* ./runregression32
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
	MATRIX x; //data set
	VECTOR y; //label set
	MATRIX xast; //data set convertito
	int n; //numero di punti del data set
	int d; //numero di dimensioni del data set
	int k; //dimensione del batch
	int degree; //grado del polinomio
	type eta; //learning rate
	//STRUTTURE OUTPUT
	VECTOR theta; //vettore dei parametri
	int t; //numero di parametri, dimensione del vettore theta
	int iter; //numero di iterazioni
	int adagrad; //accelerazione adagrad
	int silent; //silenzioso
	int display; //stampa risultati
} params;

/*
*
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere
* 	memorizzate mediante array di array (double**).
*
* 	In entrambi i casi il candidato dovr?? inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente ?? che le matrici siano in row-major order.
*
*/

void* get_block(int size, int elements) {
	return _mm_malloc(elements*size,16);
}

void free_block(void* p) {
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

/*
*
* 	load_data
* 	=========
*
*	  Legge da file una matrice di N righe						
* 	e M colonne e la memorizza in un array lineare in row-major order
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri doubleing-point a precisione singola
*
*****************************************************************************
*	Se lo si ritiene opportuno, ?? possibile cambiare la codifica in memoria
* 	della matrice.
*****************************************************************************
*
*/

MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;

	fp = fopen(filename, "rb");

	if (fp == NULL){
	printf("'%s': bad data file name!\n", filename);
	exit(0);
	}

	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}


/*
* 	save_data
* 	=========
*
*	  Salva su file un array lineare in row-major order
*	  come matrice di N righe e M colonne
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o doubleing-point a precisione singola
*/

void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	fclose(fp);
}

// PROCEDURE ASSEMBLY
int fattoriale (int k) {
	int f=1;
	for (int c = 1; c <= k; c++)
		f = f * c;
	return f;
}

int dim(int n, int k) {
	int x=0;
	int i;
	for(i=1;i<=k;i++) {
		int numeratore=fattoriale(n+i-1);
		int denominatore=fattoriale(i)*fattoriale(n-1);
		x+=numeratore/denominatore;
	}
	return x;
}

int dimNow(int n, int k) {
	int x=0;
	int numeratore=fattoriale(n+k-1);
	int denominatore=fattoriale(k)*fattoriale(n-1);
	x+=numeratore/denominatore;
	return x;
}

int dimPadding(int dimensioneReale){
	int d=dimensioneReale/16;
	int r=dimensioneReale%16;
	if(r==0) return dimensioneReale;
	else return (d+1)*16;
}

void stampaOsservazioni(double*x, double* y, int n, int d){
	for(int i=0;i<n;i++){
		printf("%s\t","[");
		for(int j=0;j<d+1;j++){
			if(j!=d)
				printf("%f\t",x[i*d+j]);
			else{
				printf("%s\t","]");
				printf("%s","[");
				printf("%f",y[i]);
				printf("%s\n","]");
			}
		}
	}
}

void stampaMatrice(int* mat, int deg, int d){
	printf("%s","[   ");
	printf("%i",mat[0]);
	printf("%s","   ");
	printf("%s\n","]");
	int i,tmp;
	int j=1, dimR=1;
	while(dimR<=deg){
		tmp=dimNow(d,dimR);
		while(tmp>0){
			printf("%s","[   ");
			for(i=0;i<dimR;i++){
				printf("%i",mat[j+i]);
				printf("%s","   ");
			}
			printf("%s\n","]");
			j+=dimR;
			tmp--;
		}
		dimR++;
	}
}

void stampaXSTAR(double* mat, int n, int d, int deg){
	int l=dim(d,deg)+1;
	int padd=dimPadding(l);
	for(int i=0;i<n;i++){
		printf("%s","[  ");
		for(int j=0;j<padd;j++){
			printf("%f",mat[i*padd+j]);
			printf("%s","   ");
		}
		printf("%s\n","]");
	}
}

void stampaTheta(double* theta, int d, int deg){
	int x=dimPadding(dim(d,deg)+1);
	int i;

	for(i=0;i<x;i++){
		printf("%s","[       ");
		printf("%.6f\t",theta[i]);
		printf("%s\n","       ]");
	}
}

void stampaArray(double* theta, int lenght){
	printf("%s\n","");
	int i;

	for(i=0;i<=lenght;i++){
		printf("%s","[       ");
		printf("%f\t",theta[i]);
		printf("%s\n","       ]");
	}
}

void stampaEQM(params* input, double* theta){
	int lenght=dim(input->d,input->degree)+1;
	int padd=dimPadding(lenght);
	double* y=input->y;
	double* xStar=input->xast;
	int n=input->n;
	double diff=0;
	for(int i=0; i<n;i++){
		double f=0.0;
		for(int j=0;j<=lenght;j++)
			f+=theta[j]*xStar[i*padd+j];
		f=y[i]-f;
		diff+=f*f;
	}
	
	printf("%s\n","L'errore quadratico medio e': ");
	printf("%f\n",diff/n);
}

int stampaLunghezza(int deg, int d){
	int lunghezza=1;
	int i,tmp;
	int j=1, dimR=1;
	while(dimR<=deg){
		tmp=dimNow(d,dimR);
		lunghezza+=tmp*dimR;
		dimR++;
	}
	return lunghezza;
}

int indice=1;
void xi_star(int * x, int * pos, int num_corr, int h, int at, int d, int* mat){
	//popolo matrice 
	int i;

	if (num_corr == h) { //condizione di uscita della ricorsione
		if (!pos) return;

		for (i = 0; i < h; i++){
			mat[indice+i]=x[pos[i]];
		}
	indice+=h;
	return;
	}

	//at ? l'indice che ricorda l'ultimo valore da cui iniziare per la prossima combinazione
	//num_corr indica l'iterazione attuale
	
	for (i = at; i < d; i++) {
		if (pos) // MECCANISMO PER NON ANDARE IN CORE DUMP
		pos[num_corr] = i;
		xi_star(x,pos, num_corr + 1, h, i, d,mat);
	}	
}

int* genera_xStar(int d, int degree){
	//riempie l'intera matrice allocata con gli indici delle combinazioni per x*
	int* coordinate=(int*)get_block(sizeof(int), d);
	int lunghezza=dim(degree,d)+1;
	int h;
	for(h=1;h<=d;h++){
		coordinate[h-1]=h; // genera [1, 2, ...., n] che corrispondono a x1, x2, ..., xn
	}
	int pos[degree]; // pari al massimo valore che assumer?? quindi pari a degree
	int* m=(int*)get_block(sizeof(int),stampaLunghezza(degree,d));
	m[0]=0;
	//richiamo x_star per ogni blocco di dimNow elementi da riempire
	
	for(h=1;h <= degree; h++){
		xi_star(coordinate,pos, 0, h, 0, d,m);
	}
	return m;
}

void calcolaRigaXStar(int* mat, double* arr, int posizione, int degree, int d, double* xStarRes, int lenght){
	/*    
	Utilizzo la matrice degli indici ricevuta per sotuire i valori
	di una specifica osservazione
	Ad esempio data la riga delle osservazioni:                  x1		     x2
										[	2.000000	2.000000	]
	
					0	        x1	         x2 		x1*x1	x1*x2	  x2*x2
	Ottengo la riga: [  1.000000   2.000000   2.000000   4.000000   4.000000   4.000000   ]
	*/
	int dimensioneCorrente=1;
	int j=1,indice=1;														//Aggiunta indice per vettorizzazione
	int tmp,i;
	double prod;

	int pos=posizione*lenght;
	int pos2=posizione*d-1;
	xStarRes[pos]=1.0;
	for (dimensioneCorrente; dimensioneCorrente<=degree; dimensioneCorrente++) {
	/*Calcola quanti elementi ha la dimensione corrente
	Ovvero quante righe scorrere nella matrice degli indici x*	  
	*/
		tmp=dimNow(d,dimensioneCorrente);
		while(tmp>0){
			prod=1.0;
			/*Il ciclo for cicla gli elementi di una singola riga della matrice degli indici
			Ovvero una singola conbinazioni degli indici*/
			for(i=0;i<dimensioneCorrente;i++){
				prod=prod*arr[pos2+mat[indice+i]];
			}
			xStarRes[pos+j]=prod;
			j++;
			indice+=dimensioneCorrente;
			tmp--;
		}
	}
	while(j<lenght){
		xStarRes[pos+j]=0.0;
		j++;
	}	
}

double* calcolaValoriXStar(int* indici, double*x, int n, int d, int degree){
	int lenght = dim(d,degree)+1;
	int lunghezza=dimPadding(lenght);
	double* xStarRes=(double*)get_block(sizeof(double),n*lunghezza);
	int i;
	/*
	Abbiamo memorizzato tutte le combinazioni di coefficienti ( x* ) dentro
	la matrice xStarRes. Ora xStarRes[i] "linearizza" la matrice per poi compiere
	una sostituzione per ognuna delle entry del dataset
	Ad esempio se xStar[i]= [1, 1] (coefficente x1*x1), con x1=5 l'osservazione 1 della riga 
	i-esima. Calcolo e inserisco nella nuova matrice 5*5
	*/
	 /*Ogni iterazione del for calcola xi* ovvero sostituisce ad x* i valori
         dell'array x corrente */
	
	for(i=0;i<n;i+=4){
		calcolaRigaXStar(indici, x, i, degree, d, xStarRes, lunghezza);
		calcolaRigaXStar(indici, x, i+1, degree, d, xStarRes, lunghezza);		//meccanismo di LOOP Unrolling
		calcolaRigaXStar(indici, x, i+2, degree, d, xStarRes, lunghezza);
		calcolaRigaXStar(indici, x, i+3, degree, d, xStarRes, lunghezza);
	}
	
	i=i-4;
	
	//spingo la macchina ad eseguire operazioni contigue in parallelo
	//dato che non ci sono dipendenze tra i valori
	
	for(;i<n;i++){
		calcolaRigaXStar(indici, x, i, degree, d, xStarRes, lunghezza);
	}
	return xStarRes;
}

double* convert_data(params* input){
	// -------------------------------------------------
	// Codificare qui l'algoritmo di conversione
	// -------------------------------------------------
	int* m = genera_xStar(input->d, input->degree);
	return calcolaValoriXStar(m, input->x, input->n, input->d, input->degree);
}

extern void batch64Prod(int p, int pfin, int i, double* osservazioni, int lenght, double* y, double* tmp, double* theta); 

void BmetodoOPM2(int p, double* theta, double* osservazioni, int lenght, double* tmp, double* y, int i){
	//batch64Prod(p,p+2,i, osservazioni, lenght,y,tmp,theta);
	for(int x = p; x<p+2;x++){
		double prodScal=0;
		//prodotto scalare
		for(int j=0;j<lenght;j++){
			prodScal+=theta[j]*osservazioni[x*lenght+j];
		}
		//sottraggo yi
		prodScal=prodScal-y[x];

		//moltiplico per xi*
		for(int j=0;j<lenght;j++){
			tmp[(x-i)*lenght+j]=osservazioni[x*lenght+j]*prodScal;
		}
	}
}

void BmetodoOPM1(int p, double* theta, double* osservazioni, int lenght, double* tmp, double* y, int i){
	//batch64Prod(p,p+1,i, osservazioni, lenght, y, tmp, theta); 
	double prodScal=0;		
	int x=p;
	for(int j=0;j<lenght;j++){
		prodScal+=theta[j]*osservazioni[x*lenght+j];
	}
	//sottraggo yi
	prodScal=prodScal-y[x];
	
	//moltiplico per xi*
	for(int j=0;j<lenght;j++){
		indice = (x-i)*lenght+j;
		tmp[indice]=osservazioni[x*lenght+j]*prodScal;
	}
}

void Bmetodo2(int i, int v, double* theta, double*osservazioni, int lenght, double*tmp, double*y){
	#pragma omp paraller for
	for(int p=i;p<=i+v-2;p+=2){
		BmetodoOPM2(p, theta, osservazioni, lenght, tmp, y, i);
	}		
}

void Bmetodo1(int indice, int i, int v, double* theta, double*osservazioni, int lenght, double*tmp, double*y){
	BmetodoOPM1(indice, theta, osservazioni, lenght, tmp, y, i);
}

double* sgdBatch(params* input, int lenght){

	//start sgd
	double* osservazioni=input->xast;
	lenght=dimPadding(lenght);
	double* theta=(double*)get_block(sizeof(double), lenght);
	int it=0;
	int i;
	double* y=input->y;
	int iter=input->iter;
	int n=input->n;
	double rate=input->eta;
	int k=input->k;
	
	double* tmp=(double*)get_block(sizeof(double),k*lenght);

	for(int i=0;i<lenght;i++){
		theta[i]=0.0;
	}
	while(it<iter){
		
		for(i=0; i<n; i+=k){ //FOR dei passi di batch

			for(int l=0;l<k*lenght;l++)
				tmp[l]=0.0;

			int j;
			int v;
			if(n-i>k)
				v=k;
			else
				v=n-i;
			indice=i+v-1;
			
			Bmetodo2(i, v, theta, osservazioni, lenght, tmp, y);
			Bmetodo1(indice, i, v, theta, osservazioni, lenght, tmp, y);
			
			for(int j=0;j<lenght;j++){
				for(int i=1;i<v;i++){
					tmp[j]+=tmp[i*lenght+j];
				}
			}

			double rapporto = rate/v;

			//sottraggo a theta
			for(j=0;j<lenght;j++){
				theta[j]=theta[j]-tmp[j]*rapporto;
			}
		}
		it++;
	}
	return theta;
}

extern void adagrad64Prod(int p, int pfin, double* theta, double* osservazioni, int lenght, double* gj, double* Gj, double* y, int i);
extern void aggiornamentoG(int p, int pfin, int lenght, double* sommatoria, double rate, double* Gj, double* gj);

void metodoOPM2(int p, double* theta, double* osservazioni, int lenght, double* gj, double* Gj, double* y, int i){
	//adagrad64Prod(p,p+2,theta, osservazioni, lenght, gj, Gj,y,i); 
	for(int x = p; x<p+2;x++){
		double prodScal=0;
		//prodotto scalare
		for(int j=0;j<lenght;j++){
			prodScal+=theta[j]*osservazioni[x*lenght+j];
		}
		//sottraggo yi
		prodScal=prodScal-y[x];

		//moltiplico per xi*
		for(int j=0;j<lenght;j++){
			indice=(x-i)*lenght+j;
			gj[indice]=osservazioni[x*lenght+j]*prodScal;
			Gj[indice]+=gj[indice]*gj[indice];
		}
	}
}

void metodoOPM1(int p, double* theta, double* osservazioni, int lenght, double* gj, double* Gj, double* y, int i){
	//adagrad64Prod(p,p+1,theta, osservazioni, lenght, gj, Gj,y,i); 
	double prodScal=0;
	//prodotto scalare
	for(int j=0;j<lenght;j++){
		prodScal+=theta[j]*osservazioni[p*lenght+j];
	}
	//sottraggo yi
	prodScal=prodScal-y[p];

	//moltiplico per xi*
	for(int j=0;j<lenght;j++){
		indice=(p-i)*lenght+j;
		gj[indice]=osservazioni[p*lenght+j]*prodScal;
		Gj[indice]+=gj[indice]*gj[indice];
	}
}

void metodo2(int i, int v, double* theta, double*osservazioni, int lenght, double*gj, double* Gj, double*y){
	#pragma omp paraller for
	for(int p=i;p<=i+v-2;p+=2){
		metodoOPM2(p, theta, osservazioni, lenght, gj, Gj, y, i);
	}		
}

void metodo1(int indice, int i, double* theta, double*osservazioni, int lenght, double*gj, double* Gj, double*y){
	metodoOPM1(indice, theta, osservazioni, lenght, gj, Gj, y, i);
}

double* sgdAdagrad(params* input, int lenght){

	//start sgd
	double* osservazioni=input->xast;
	lenght=dimPadding(lenght);
	double* theta=(double*)get_block(sizeof(double), lenght+100);
	int i;
	double* y=input->y;
	int iter=input->iter;
	int n=input->n;
	double rate=input->eta;
	int k=input->k;
	double eps=1E-8;
	int indice;

	double* Gj=(double*)get_block(sizeof(double),k*lenght);
	double* gj=(double*)get_block(sizeof(double),k*lenght);
	double* sommatoria=(double*)get_block(sizeof(double),lenght);

	for(int p1=0;p1<k;p1++)
		for(int p2=0;p2<lenght;p2++){
			Gj[p1*lenght+p2]=0.0;
			gj[p1*lenght+p2]=0.0;
		}

	for(int i=0;i<lenght;i++)
		theta[i]=0.0;

	for(int it=0;it<iter;it++){
		for(i=0; i<n; i=i+k){

			int j;
			int v;
			if(n-i>k)
				v=k;
			else
				v=n-i;
			
			indice=i+v-1;
			
			metodo2(i, v, theta, osservazioni, lenght, gj, Gj, y);
			metodo1(indice, i, theta, osservazioni, lenght, gj, Gj, y);
			
			for(int j=0;j<lenght;j++)
				sommatoria[j]=0.0;
			
			int p=i;
			int pfin=i+v-1;
			
			for(int p=i; p<i+v;p++)
				for(int j=0;j<lenght;j++){
					indice = ((p-i) * lenght) +j;
					sommatoria[j]+=(rate/sqrt(Gj[indice]+eps))*gj[indice];
				}
			
			
			//sottraggo a theta
			for(j=0;j<lenght;j++){
				theta[j]=theta[j]-sommatoria[j]/v;
			}
		
		}

	}
	return theta;
}

int main(int argc, char** argv) {

	double* osservazioni;	
	char fname[256];
	char* dsname;
	char* filename;
	int i, j, k;
	clock_t t;
	double time;
	int yd = 1;

	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->x = NULL;
	input->y = NULL;
	input->xast = NULL;
	input->n = 0;
	input->d = 0;
	input->k = -1;
	input->degree = -1;
	input->eta = -1;
	input->iter = -1;
	input->adagrad = 0;
	input->theta = NULL;
	input->t = 0;
	input->adagrad = 0;
	input->silent = 0;
	input->display = 0;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s D -batch <k> -degree <deg> -eta <eta> -iter <it> [-adagrad]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tD: il nome del file, estensione .data per i dati x, estensione .labels per le etichette y\n");
		printf("\t-batch <k>: il numero di campini nel batch\n");
		printf("\t-degree <deg>: il grado del polinomio\n");
		printf("\t-eta <eta>: il learning rate\n");
		printf("\t-iter <it>: il numero di iterazioni\n");
		printf("\t-adagrad: l'acceleratore AdaGrad\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (par == 1) {
			filename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-batch") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing batch dimension value!\n");
				exit(1);
			}
			input->k = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-degree") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing degree value!\n");
				exit(1);
			}
			input->degree = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-eta") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing eta value!\n");
				exit(1);
			}
			input->eta = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-iter") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing iter value!\n");
				exit(1);
			}
			input->iter = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-adagrad") == 0) {
			input->adagrad = 1;
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if(filename == NULL || strlen(filename) == 0){
		printf("Missing input file name!\n");
		exit(1);
	}

	dsname = basename(strdup(filename));
	sprintf(fname, "%s.data", filename);
	input->x = load_data(fname, &input->n, &input->d);
	sprintf(fname, "%s.labels", filename);
	input->y = load_data(fname, &input->n, &yd);
	sprintf(fname, "%s.theta", filename);
	
	if(input->k < 0){
		printf("Invalid value of batch dimension parameter!\n");
		exit(1);
	}

	if(input->degree < 0){
		printf("Invalid value of degree parameter!\n");
		exit(1);
	}

	if(input->eta < 0){
		printf("Invalid value of eta parameter!\n");
		exit(1);
	}

	if(input->iter < 0){
		printf("Invalid value of iter parameter!\n");
		exit(1);
	}

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Input data name: '%s.data'\n", filename);
		printf("Input label name: '%s.labels'\n", filename);
		printf("Data set size [n]: %d\n", input->n);
		printf("Number of dimensions [d]: %d\n", input->d);
		printf("Batch dimension: %d\n", input->k);
		printf("Degree: %d\n", input->degree);
		printf("Eta: %f\n", input->eta);
		if(input->adagrad)
			printf("Adagrad enabled\n");
		else
			printf("Adagrad disabled\n");
	}
	
	//
	// Conversione Dati
	//

	t = clock();
	input->xast=convert_data(input);
	t = clock() - t;
	time = ((double)t)/CLOCKS_PER_SEC;
	sprintf(fname, "%s.xast", dsname);

	if(!input->silent)
		printf("Conversion time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	int dimensioneTheta=dim(input->d, input->degree)+1;
	double* theta;
	
	//
	// Regressione
	//

	t = clock();
	if(!input->adagrad){
		sprintf(fname, "%s.theta.sgdomp", dsname);
		theta=sgdBatch(input, dimensioneTheta);
	}
	else{
		sprintf(fname, "%s.theta.adagradomp", dsname);
		theta=sgdAdagrad(input, dimensioneTheta);
	}
	t = clock() - t;
	time = ((double)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("Regression time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato di theta
	//
	
	save_data(fname, theta, dimensioneTheta, 1);
	
	if(input->display){
		printf("theta: [");
		for(i=0; i<input->t-1; i++)
		printf("%f,", input->theta[i]);
		printf("%f]\n", input->theta[i]);
	}
	
	if(!input->silent){
		if(!input->adagrad){
			printf("%s\n","				Vettore Theta sgdBatch");
			stampaTheta(theta, input->d,input->degree);
			stampaEQM(input, theta);
		}else{
			printf("%s\n","				Vettore Theta sgdAdagrad");
			stampaTheta(theta, input->d,input->degree);
			stampaEQM(input, theta);
		}
		printf("\nDone.\n");
	}
	return 0;
}
