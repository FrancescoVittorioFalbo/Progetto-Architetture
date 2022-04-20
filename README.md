# Progetto-Architetture
Progetto di gruppo per il corso di Architetture e Programmazione dei Sistemi di Elaborazione

---


**Caso specifico**
“Polynomial Regression and Stochastic Gradient Descent”
in linguaggio assembly x86-32+SSE, x86-64+AVX e openMP
**Descrizione problema**

Dato un insieme di osservazioni concernenti le variabili x e y e fissata una famiglia di funzioni
f(x, Θ) da usare per approssimare la relazione tra x e y, ossia
y ∼ f(x, Θ)
il problema consiste nel trovare i parametri Θ della funzione f per cui sia minimizzato l’errore
di approssimazione. Θ = [θ1, . . . , θt] `e in generale un vettore di t parametri
In particolare, in questo progetto assumiamo che la funzione f sia una funzione polinomiale di
grado deg e che l’algoritmo di ricerca dei parametri si basi sulla discesa stocastica del gradiente.
