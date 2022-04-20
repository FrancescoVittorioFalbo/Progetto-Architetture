# Progetto-Architetture
Progetto realizzato in gruppo per il corso "Architetture e Programmazione dei Sistemi di Elaborazione"

---

Obiettivo: Sviluppare in linguaggio assembly x86-32+SSE e x86-64+AVX e openMP l'algoritmo "Polynomial Regression and Stochastic Gradient Descent"

---

Dato un insieme di osservazioni concernenti le variabili _x_ e _y_ e fissata una famiglia di funzioni _f(x,Θ)_ da usare per approssimare la relazione tra _x_ e _y_, ossia <div align="center"> <i>y∼f(x,Θ)</i></div>

il problema consiste nel trovare i parametri _Θ_ della funzione _f_ per cui sia minimizzato l’errore di approssimazione. _Θ_ = [_θ<sub>1</sub>_, ..., _θ<sub>t</sub>_] è in generale un vettore di _t_ parametri.
In particolare, in questo progetto assumiamo che la funzione _f_ sia una funzione polinomiale di grado _deg_ e che l’algoritmo di ricerca dei parametri si basi sulla discesa stocastica del gradiente.


**Preliminari**

Un’osservazione è una coppia (_x_, _y_) dove _x_ = [_x<sub>1</sub>_, ..., _x<sub>d</sub>_] ∈ _R_<sup>d</sup> è un vettore _d_ dimensionale e _y_ ∈ _R_ è un reale.
Dato un insieme di osservazioni _D_ = {(_x<sub>1</sub>_, _y<sub>1</sub>_),... ,(_x<sub>n</sub>_, _y<sub>n</sub>_)}, |_D_| = _n_ indica la cardinalità dell’insieme _D_, ossia il numero di osservazioni presenti in _D_.

**Regressione Polinomiale**

Dato un insieme di osservazioni concernente la variabile _d_ dimensionale _x_ ∈ R<sup>d</sup> e la variabile uni-dimensionale _y_ ∈ R, la _regressione_ rappresenta un metodo di stima della variabile dipendente _y_ dato il valore della variabile multi dimensionale indipendente _x_ = [_x<sub>1</sub>_, ..., _x<sub>d</sub>_] nell’ipotesi che _f_ sia una funzione polinomiale di un certo grado _deg_, ovvero:

<p align="center">
  <img src=https://i.postimg.cc/k5GQTZvF/1.png" />
</p>

dove [_d_] rappresenta l’insieme {1, ..., _d_} e [_d_]<sup>h</sup> indica il prodotto cartesiano

<p align="center">
  <img src="https://i.postimg.cc/dtDJtj6N/2.png" />
</p>

L’equazione (1) può essere rappresentata in forma vettoriale

<p align="center">
  <img src="https://i.postimg.cc/MHQMBm5w/3.png" />
</p>

dove ⟨·, ·⟩ indica il prodotto scalare e, dato un vettore **x** _d_ dimensionale, **x**<sup>\*</sup> indica il vettore i cui elementi **x**<sup>\*</sup><sub>j</sub> sono ottenuti come

<p align="center">
  <img src="https://i.postimg.cc/P5vh1bC3/4.png" />
</p>

Quindi, ad esempio, nel caso in cui **x** = [_x_] sia una variabile uni-dimensionale, ossia _d_ = 1, e _deg_ = 4 si ottiene

<p align="center">
  <img src="https://i.postimg.cc/J4ZCnrvv/5.png" />
</p>

o, rinominando i parametri,

<p align="center">
  <img src="https://i.postimg.cc/pVCyXt52/6.png" />
</p>

o in forma vettoriale

<p align="center">
  <img src="https://i.postimg.cc/P5pq4p3x/7.png" />
</p>

Nel caso in cui _x_ = [_x_<sub>1</sub>, _x_<sub>2</sub>] sia una variabile bi-dimensionale, ossia _d_ = 2, e _deg_ = 4 si ottiene

<p align="center">
  <img src="https://i.postimg.cc/nzZcpdHj/8.png" />
</p>

o, rinominando i parametri,

<p align="center">
  <img src="https://i.postimg.cc/pX3KXNy1/9.png.png" />
</p>

o in forma vettoriale

<p align="center">
  <img src="https://i.postimg.cc/Sxpjvmnq/10.png" />
</p>

con

<p align="center">
  <img src="https://i.postimg.cc/3xMDQgbs/11.png" />
</p>

e

<p align="center">
  <img src="https://i.postimg.cc/y8F31kMk/12.png" />
</p>

Si noti che ci sono termini omessi, perchè generano termini accorpabili a quelli espressi, si consideri ad esempio il caso _h_ = 3, si avrebbe

<p align="center">
  <img src="https://i.postimg.cc/yNjJR2YK/13.png" />
</p>

e quindi

<p align="center">
  <img src="https://i.postimg.cc/43fb7DRv/14.png" />
</p>

semplificabile come

<p align="center">
  <img src="https://i.postimg.cc/nL3Kv3ny/15.png" />
</p>
e θ<sub>1,1,2</sub> + θ<sub>2,1,1</sub> + θ<sub>1,2,1</sub> è di fatto un unico parametro così come θ<sub>1,2,2</sub> + θ<sub>2,1,2</sub> + θ<sub>2,2,1</sub>.

Il problema da risolvere è quindi specificato dalla Definizione 1.

**Definizione 1** _Dato un data set D_ = [**x**_, y_] _e un grado_ deg _, trovare i parametri_ Θ _della funzione f_(**x**, Θ) _specificata dall’Equazione 1 per cui sia minimizzato l’errore quadratico medio calcolato come:_

<p align="center">
  <img src="https://i.postimg.cc/k4ZxHD6w/16.png" />
</p>

**Discesa Stocastica del Gradiente**

La discesa stocastica del gradiente è un metodo iterativo che, ad ogni iterazione, sostituisce il valore esatto del gradiente della funzione costo con una stima ottenuta valutando il gradiente solo su un sottinsieme degli addendi. È ampiamente usato per l’allenamento di una varietà di modelli di apprendimento automatico, come macchine a vettori di supporto, regressione logistica e modelli grafici. Più in detttaglio, l’algoritmo stima i parametri Θ che minimizzano una funzione di costo _C_(Θ, **x**, _y_), aggiornando iterativamente Θ e, in particolare, Θ<sub>t+1</sub> viene ottenuto sottraendo a Θ<sub>t</sub> il gradiente della funzione di costo _C_ calcolato rispetto a Θ, dove Θ<sub>t</sub> indica il valore dei parametri all’interazione _t_. In formula:

<p align="center">
  <img src="https://i.postimg.cc/MpgW1hsj/17.png" />
</p>

Nel caso particolare di funzione di costo pari all’errore quadratico medio,

<p align="center">
  <img src="https://i.postimg.cc/7L7q9mct/18.png" />
</p>

Lo pseudo codice relativo all’algoritmo, considerando regressione polinomiale ed errore quadratico medio come funzione di costo, è descritto dalla funzione SGD, la variante utilizzante un batch di _k_ campioni è descritto dalla funzione SGDBATCH e la variante utilizzante l’algoritmo _AdaGrad_ come acceleratore è descritto dalla funzione AdaGrad.

<p align="center">
  <img src="https://i.postimg.cc/FRVN51hj/19.png" />
</p>

#
# Descrizione dell’attività progettuale
Obiettivo del progetto è mettere a punto un’implementazione dell’algoritmo di Regressione Polinomiale in linguaggio C e di migliorarne le prestazioni utilizzando le tecniche di ottimizzazione basate sull’organizzazione ell’hardware.
L’ambiente sw/hw di riferimento è costituito dal linguaggio di programmazione C (gcc), dal linguaggio assembly x86-32+SSE e dalla sua estensione x86-32+AVX (nasm) e dal sistema operativo Linux (ubuntu).

<p align="center">
  <img src="https://i.postimg.cc/wMWzMPyM/20.png" />
</p>

In particolare il codice deve consentire di trovare i parametri relativi alla regressione polinomiale in un insieme di osservazioni dato in input con l’algoritmo SGDBATCH (parametro -batch <k>), con eventuale accelerazione AdaGrad, dato il grado del polinomio (parametro -degree <deg>), il valore di η (parametro -eta <eta>) e la specifica _adagrad_ (parametro -adagrad).
Quindi la chiamata avrà la seguente struttura:

<p align="center">
  <img src="https://i.postimg.cc/HnWHpHLF/21.png" />
</p>

si noti che il parametro **adagrad** è opzionale e la sua presenza indica che l’algoritmo deve prevedere l’accelerazione AdaGrad. Qualora un valore di un parametro (sia esso di default o specificato dall’utente) non sia applicabile, il codice deve segnalarlo con un messaggio e terminare.

<p align="center">
  <img src="https://i.postimg.cc/gkJdQYQZ/22.png" />
</p>

* Sono richieste due soluzioni software, la prima per l’architettura x86-32 SSE e la seconda per l’architettura x86-64+AVX.

* È inoltre richiesta per ogni soluzione, una versione che faccia uso delle istruzioni **OpenMP**. I nomi dei relativi file dovranno contenere il suffisso “\_omp” (es. regression32 omp.c).
