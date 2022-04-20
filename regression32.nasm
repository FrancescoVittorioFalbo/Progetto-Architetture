%include "sseutils32.nasm"

section .data
	dim 		equ 4
	unroll 	equ 4
	dimB 		equ 4
	elementiunrollingB	equ	64
	elementiunrolling 	equ	64
	align 16
	epsilon  	dd 0.00000001,0.00000001,0.00000001,0.00000001
	
section .bss
	;------------------------------------------->Variabili Adagrad
	alignb 16
	k		resd 1
	alignb 16
	n		resd 1
	alignb 16
	lenght	resd 1	;numero di coordinate theta per ogni oggetto
	alignb 16
	rate 		resd	1
	alignb 16
	rapporto	resd 4
	alignb 16
	indice	resd 1
	
	;----------------------------------------------
	; 	VARIABILI AUSILIARIE
	;----------------------------------------------
	
	alignb 16
	v		resd 4
	alignb 16
	limite 	resd 1	;contiene i+v
	alignb 16
	i		resd	1
	alignb 16		
	p		resd	1
	;------------------------------------------------>Variabili Batch
	alignb 16
	kB		resd 1
	alignb 16
	nB		resd 1
	alignb 16
	lengthB	resd 1	;numero di coordinate theta per ogni oggetto
	alignb 16
	rateB 	resd	1
	alignb 16
	rapportoB	resd 4
	alignb 16
	tmpB	resd 1
	
	;----------------------------------------------
	; 	VARIABILI AUSILIARIE
	;----------------------------------------------
	
	alignb 16
	vB		resd 1
	alignb 16
	auxB		resd 1
	alignb 16
	limiteB 	resd 1
	alignb 16
	iB		resd	1
	alignb 16		
	pB		resd	1
section .text
	
	input			equ 8
	osserv			equ 12
	theta 			equ 16
	startlenght 		equ 20
	Gj				equ 24
	starttmpB			equ 24
	gj				equ 28
	sommatoria		equ 32

global adagrad32
global batch32

adagrad32:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	ebp							; salva il Base Pointer
		mov		ebp, esp						; il Base Pointer punta al Record di Attivazione corrente
		push	ebx							; salva i registri da preservare
		push	esi
		push	edi							
		
		
		MOV	EAX, [EBP + startlenght]
		MOV	[lenght], EAX
		
		MOV	EAX, [EBP+input]
		MOV	EBX, [EAX + 12]
		MOV	[n], EBX
		MOV	EBX, [EAX + 20]
		MOV	[k], EBX
		MOV	EBX, [EAX + 28]
		MOV	[rate], EBX
		
		MOV 	EBX, 0					;EBX indice i del for dei passi di batch
fori:
		
		MOV	[i], EBX
		MOV 	ECX, [n]		 ;In ECX inserisco il valore v, contatore di elementi nel corrente passo di batch
		SUB 	ECX, EBX
		CMP 	ECX, [k]
		
		JLE		salto
		
		MOV 	ECX, [k]
salto:
		CVTSI2SS	XMM0, ECX		;convertire integer to float
		SHUFPS	XMM0, XMM0, 00000000
		MOVAPS	[v], XMM0
		
		MOV 	EDX, EBX ;EDX indice p per scorrere gli elementi del corrente passo di batch faccio p=i
		MOV 	ESI, EBX
		ADD 	ESI, ECX ;In ESI ho inserito i+v per il forp
		;edx ho p e in esi ho i+v condizione di inizio del forp
		
		MOV	[limite], ESI
		
forp:
		MOV 	[p], EDX
		
		XORPS	XMM1, XMM1				;XMM1 contiene prodScal
		MOV 	EDI, 0					;EDI indice j del for del prodotto scalare
		MOV 	ESI, 0
		MOV 	ECX, EDX				;metto in ECX il valore di p
		IMUL	ECX, [lenght]				;moltiplico per lenght
		IMUL	ECX, ECX, dim				;IMUL ci salva p*lenght*dim
		
		
		MOV	EBX,[EBP + theta]
		MOV	EDX,[EBP + osserv]
		
forj1:	
		MOVAPS	XMM2, [EBX + EDI]
		MULPS	XMM2, [EDX + ECX]
		ADDPS	XMM1, XMM2
		MOVAPS	XMM3, [EBX + EDI + 16]
		MULPS	XMM3, [EDX + ECX + 16]
		ADDPS	XMM1, XMM3
		MOVAPS	XMM4, [EBX + EDI + 32]
		MULPS	XMM4, [EDX + ECX + 32]
		ADDPS	XMM1, XMM4
		MOVAPS	XMM5, [EBX + EDI + 48]
		MULPS	XMM5, [EDX + ECX + 48]
		ADDPS	XMM1, XMM5	
		
		;questo chiude il forj1
		ADD		ESI, 16 
		ADD		EDI, elementiunrolling
		ADD 	ECX, elementiunrolling
		CMP	ESI, [lenght] 
		JL		forj1

		HADDPS	XMM1, XMM1
		HADDPS	XMM1, XMM1				;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		MOV	EDX, [EBP + input]
		MOV	EBX, [EDX+4]
		
		MOV 	EDX, [p]
		
		IMUL	EDI, EDX, dim
		SUBSS	XMM1, [EBX + EDI]		;in XMM1 ho prodScal = prodScal - y[p]
		SHUFPS	XMM1, XMM1, 00000000
		
		MOV 	EDI, 0					;EDI indice j del for del prodotto scalare
		MOV 	ESI, EDX
		SUB		ESI, [i]
		IMUL	ESI, [lenght]
		IMUL 	ESI, ESI, dim				;calcolo indice = p-i * lenght * dim
		
		
		IMUL	EDX, [lenght]				;posso riutilizzare EDX tanto il valore dell'indice p e salvato nella variabile p
		IMUL	EDX, EDX, dim			;IMUL ci salva p*lenght*dim in EDX

		
		MOV 	EAX, [EBP + gj]
		MOV	EBX, [EBP + osserv]
		MOV	ECX, [EBP + Gj]
		
forj2:	
		MOVAPS	XMM0, XMM1				;in aux ce prodScal
		MOVAPS	XMM2, XMM1
		MOVAPS	XMM4, XMM1
		MOVAPS	XMM6, XMM1
			
		MULPS	XMM0, [EBX + EDX]
		MULPS	XMM2, [EBX + EDX + 16]
		MULPS	XMM4, [EBX + EDX + 32]
		MULPS	XMM6, [EBX + EDX + 48]		
		
		MOVAPS	[EAX + ESI], XMM0		;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		MOVAPS	[EAX + ESI + 16], XMM2
		MOVAPS	[EAX + ESI + 32], XMM4
		MOVAPS	[EAX + ESI + 48], XMM6
		
		MULPS	XMM0, XMM0
		MULPS	XMM2, XMM2
		MULPS	XMM4, XMM4
		MULPS	XMM6, XMM6
		
		ADDPS	XMM0, [ECX + ESI]
		ADDPS	XMM2, [ECX + ESI + 16]
		ADDPS	XMM4, [ECX + ESI + 32]
		ADDPS	XMM6, [ECX + ESI + 48]
		
		MOVAPS	[ECX + ESI], XMM0		;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		MOVAPS	[ECX + ESI + 16], XMM2
		MOVAPS	[ECX + ESI + 32], XMM4
		MOVAPS	[ECX + ESI + 48], XMM6
		
		;questo chiude il forj2
		ADD		EDI, 16
		ADD		ESI, elementiunrolling
		ADD		EDX, elementiunrolling
		CMP	EDI, [lenght]
		JL		forj2
		
		; questo chiude il forp
		MOV	EDX, [p]
		MOV 	ESI, [limite]
		ADD		EDX, 1
		CMP 	EDX, ESI
		JL		forp

		MOV 	ESI, 0					;ESI indirizzo di partenza dei 4 float di interesse
		MOV 	EDI, 0					;EDI numero di valori che scorriamo in un ciclo, (indice l del for di azzeramento di tmp)
		MOV	EAX, [EBP + sommatoria]
forl:
		
		XORPS 	XMM0, XMM0				;XMM0 ho sommatoria 
		MOVAPS	[ EAX + ESI ], XMM0
		MOVAPS	[ EAX + ESI  + 16], XMM0
		MOVAPS	[ EAX + ESI  + 32], XMM0
		MOVAPS	[ EAX + ESI  + 48], XMM0
		ADD 	ESI, elementiunrolling			
		ADD 	EDI, 16
		CMP 	EDI, [lenght]
		JL 		forl
		
		MOV	ESI, [i]
		MOV	EDI, [limite]
		MOV	[p], ESI
		
forp2:
		MOV 	EAX, [p]
		SUB		EAX, [i]
		IMUL	EAX, [lenght]
		IMUL 	EAX, EAX, dim				;calcolo indice = p-i * lenght * dim e lo metto in EAX
		MOV	[indice], EAX
		MOV	EBX, 0
		MOV	ECX, [EBP + gj]
		
forj3:
		MOV	[indice], EAX
		MOV	EDX, [EBP + Gj]
		IMUL	EBX, EBX, dim				; in EBX ho j
		
		MOVAPS	XMM4, [EDX + EAX]			;qui abbiamo Gj
		MOVAPS	XMM5, [EDX + EAX + 16]
		MOVAPS 	XMM6, [EDX + EAX + 32]
		MOVAPS	XMM7, [EDX + EAX + 48]
		
		ADDPS	XMM4, [epsilon]
		ADDPS	XMM5, [epsilon]
		ADDPS	XMM6, [epsilon]
		ADDPS	XMM7, [epsilon]
		
		SQRTPS	XMM4, XMM4
		SQRTPS	XMM5, XMM5
		SQRTPS	XMM6, XMM6
		SQRTPS	XMM7, XMM7
		
		MOVSS	XMM0, [rate]
		SHUFPS	XMM0, XMM0, 00000000
		
		MOVAPS 	XMM1, XMM0
		MOVAPS	XMM2, XMM0
		MOVAPS	XMM3, XMM0
		
		DIVPS	XMM0, XMM4
		DIVPS	XMM1, XMM5
		DIVPS	XMM2, XMM6
		DIVPS	XMM3, XMM7
		
		MOVAPS	XMM4, [ECX + EAX]			;qui abbiamo gj
		MOVAPS	XMM5, [ECX + EAX + 16]
		MOVAPS 	XMM6, [ECX + EAX + 32]
		MOVAPS	XMM7, [ECX + EAX + 48]
		
		
		MULPS	XMM0, XMM4
		MULPS	XMM1, XMM5
		MULPS	XMM2, XMM6
		MULPS	XMM3, XMM7
		
		MOV	EDX, [EBP + sommatoria]
		
		ADDPS	XMM0, [EDX + EBX]			;qui abbiamo sommatoria
		ADDPS	XMM1, [EDX + EBX + 16]
		ADDPS 	XMM2, [EDX + EBX + 32]
		ADDPS	XMM3, [EDX + EBX + 48]
		
		MOVAPS 	[EDX + EBX], XMM0
		MOVAPS	[EDX + EBX + 16], XMM1
		MOVAPS	[EDX + EBX + 32], XMM2
		MOVAPS	[EDX + EBX + 48], XMM3
		
		;questo chiude forj3
		MOV	EDX, 0
		MOV	EAX, EBX
		MOV	EBX, dim
		DIV		EBX
		MOV	EBX, EAX
		
		
		ADD		EBX, 16
		MOV	EAX, [indice]
		ADD		EAX, elementiunrolling
		CMP	EBX, [lenght]
		JL		forj3
		
		;questo chiude forp2
		ADD		ESI, 1
		MOV	[p],ESI	
		CMP	ESI, EDI
		JL		forp2
		
		
		MOV 	EDI, 0
		MOV 	ECX, 0
		MOV	EBX, [EBP + theta]
		MOV	EAX, [EBP + sommatoria]
		
forj4:	
		MOVAPS 	XMM2, [EAX + EDI]
		MOVAPS 	XMM3, [EAX + EDI + 16]
		MOVAPS 	XMM4, [EAX + EDI + 32]
		MOVAPS 	XMM5, [EAX + EDI + 48]
		
		DIVPS 	XMM2, [v]
		DIVPS 	XMM3, [v]
		DIVPS 	XMM4, [v]
		DIVPS 	XMM5, [v]
		
		MOVAPS 	XMM0, [EBX + EDI]
		MOVAPS 	XMM1, [EBX + EDI + 16]
		MOVAPS 	XMM6, [EBX + EDI + 32]
		MOVAPS 	XMM7, [EBX + EDI + 48]
		
		SUBPS 	XMM0, XMM2					;salvataggio in theta dei valori calcolati
		SUBPS 	XMM1, XMM3		
		SUBPS 	XMM6, XMM4		
		SUBPS 	XMM7, XMM5	
		
		MOVAPS	[EBX + EDI], XMM0
		MOVAPS 	[EBX + EDI + 16], XMM1
		MOVAPS 	[EBX + EDI + 32], XMM6
		MOVAPS 	[EBX + EDI + 48], XMM7
		
		;questo chiude il forj3
		ADD		ECX, 16
		ADD		EDI, elementiunrolling
		CMP	ECX, [lenght]
		JL		forj4
		
		;ripristino gli indici v e i in modo da trovarli aggiornati al prossimo for
		MOV 	EBX, [i]
		
		; questo chiude il fori
		ADD 	EBX, [k]
		CMP	EBX, [n]
		JL 		fori
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp								; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante
			
		

batch32:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	ebp							; salva il Base Pointer
		mov		ebp, esp						; il Base Pointer punta al Record di Attivazione corrente
		push	ebx							; salva i registri da preservare
		push	esi
		push	edi							
		
		
		MOV	EAX, [EBP + startlenght]
		MOV	[lengthB], EAX
		
		MOV	EAX, [EBP+input]
		MOV	EBX, [EAX + 12]
		MOV	[nB], EBX
		MOV	EBX, [EAX + 20]
		MOV	[kB], EBX
		MOV	EBX, [EAX + 28]
		MOV	[rateB], EBX
		
		MOV 	EAX, [EBP+starttmpB]
		
		MOV 	EBX, 0					;EBX indice i del for dei passi di batch
foriB:
		MOV	[iB], EBX
		MOV 	ESI, 0					;ESI indirizzo di partenza dei 4 float di interesse
		MOV 	EDI, 0					;EDI numero di valori che scorriamo in un ciclo, (indice l del for di azzeramento di tmp)
forlB:
		XORPS 	XMM0, XMM0				;XMM0 ho tmp 
		MOVAPS	[ EAX + ESI ], XMM0
		MOVAPS	[ EAX + ESI  + 16], XMM0
		MOVAPS	[ EAX + ESI  + 32], XMM0
		MOVAPS	[ EAX + ESI  + 48], XMM0
		ADD 	ESI, elementiunrollingB				; in caso sostituire con 64
		ADD 	EDI, 16
		CMP 	EDI, [lengthB]
		JL 		forlB
		;----------------------
		;CALCOLO DI V
		;---------------------
		MOV 	ECX, [nB]					;In ECX inserisco il valore v, contatore di elementi nel corrente passo di batch
		SUB 	ECX, EBX
		CMP	ECX, [kB]
		JLE 		saltoB
		
		MOV	ECX, [kB]
saltoB:									;In questo caso avremo ECX= n - i
		MOV	[vB], ECX
		MOV 	EDX, EBX 				;EDX indice p per scorrere gli elementi del corrente passo di batch faccio p=i
		MOV	ESI, EBX
		ADD 	ESI, ECX					;In ESI ho inserito i+v per il forp
		;edx ho p e in esi ho i+v condizione di inizio del forp
		
		MOV	[limiteB], ESI
		
forpB:
		XORPS	XMM1, XMM1				;XMM1 contiene prodScal
		MOV 	EDI, 0					;EDI indice j del for del prodotto scalare
		MOV 	ESI, 0
		MOV 	ECX, EDX				;metto in ECX il valore di p
		IMUL	ECX, [lengthB]				;moltiplico per lenght
		IMUL	ECX, ECX, dimB				;IMUL ci salva p*lenght*dimB
		
		
		MOV 	[pB], EDX
		MOV	EBX,[EBP + theta]
		MOV	EDX,[EBP + osserv]
		
forj1B:	
		MOVAPS	XMM2, [EBX + EDI]
		MULPS	XMM2, [EDX + ECX]
		ADDPS	XMM1, XMM2
		MOVAPS	XMM3, [EBX + EDI + 16]
		MULPS	XMM3, [EDX + ECX + 16]
		ADDPS	XMM1, XMM3
		MOVAPS	XMM4, [EBX + EDI + 32]
		MULPS	XMM4, [EDX + ECX + 32]
		ADDPS	XMM1, XMM4
		MOVAPS	XMM5, [EBX + EDI + 48]
		MULPS	XMM5, [EDX + ECX + 48]
		ADDPS	XMM1, XMM5	
		
		;questo chiude il forj1
		ADD		ESI, 16 
		ADD		EDI, elementiunrollingB
		ADD 	ECX, elementiunrollingB
		CMP	ESI, [lengthB] 
		JL		forj1B

		HADDPS	XMM1, XMM1
		HADDPS	XMM1, XMM1				;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		MOV	EDX, [EBP + input]
		MOV	EBX, [EDX+4]
		
		;ripristino l'indice p giusto prima di aggiornarlo
		MOV 	EDX, [pB]
		
		IMUL	EDI, EDX, dimB
		SUBSS	XMM1, [EBX + EDI]		;in XMM1 ho prodScal = prodScal - y[p]
		MOVSS	[auxB], XMM1				;in AUX salvo cio che ho appena calcolato
		
		
		MOV 	EDI, 0					;EDI indice j del for del prodotto scalare
		IMUL	EDX, [lengthB]				;posso riutilizzare EDX tanto il valore dell'indice p e salvato nella variabile p
		IMUL	EDX, EDX, dimB			;IMUL ci salva p*lenght*dimB

		
		MOV	EBX, [EBP + osserv]
		MOV	ECX, 0
		
forj2B:	
		XORPS	XMM1,XMM1
		XORPS	XMM3,XMM3
		XORPS	XMM5,XMM5
		XORPS	XMM7,XMM7
		
		MOVSS	XMM0, [auxB]
		MOVSS	XMM2, [auxB]
		MOVSS	XMM4, [auxB]
		MOVSS	XMM6, [auxB]
		
		SHUFPS	XMM0, XMM0, 00000000		
		MULPS	XMM0, [EBX + EDX]
		ADDPS	XMM1, XMM0		
		
		SHUFPS	XMM2, XMM2, 00000000
		MULPS	XMM2, [EBX + EDX + 16]
		ADDPS	XMM3, XMM2		
		
		SHUFPS	XMM4, XMM4, 00000000
		MULPS	XMM4, [EBX + EDX + 32]
		ADDPS	XMM5, XMM4		
		
		SHUFPS	XMM6, XMM6, 00000000
		MULPS	XMM6, [EBX + EDX + 48]
		ADDPS	XMM7, XMM6
		
		ADDPS	XMM1, [EAX + ECX]			;salvataggio in tmp dei valori calcolati, ovvero tutto ciï¿½ che sta a destra della sommatoria
		ADDPS	XMM3, [EAX + ECX + 16]
		ADDPS	XMM5, [EAX + ECX + 32]
		ADDPS	XMM7, [EAX + ECX + 48]
		
		
		MOVAPS	[EAX + ECX], XMM1			;salvataggio in tmp dei valori calcolati, ovvero tutto ciï¿½ che sta a destra della sommatoria
		MOVAPS	[EAX + ECX + 16], XMM3
		MOVAPS	[EAX + ECX + 32], XMM5
		MOVAPS	[EAX + ECX + 48], XMM7
		
		;questo chiude il forj2
		ADD		EDI, 16
		ADD		ECX, elementiunrollingB
		ADD		EDX, elementiunrollingB
		CMP	EDI, [lengthB]
		JL		forj2B
		
		; questo chiude il forp
		MOV	EDX, [pB]
		MOV 	ESI, [limiteB]
		ADD		EDX, 1
		CMP 	EDX, ESI
		JL		forpB


		MOVSS	XMM0, [rateB]	
		CVTSI2SS	XMM1, [vB]						;CAPIRE COME LEGGERE V AL POSTO DI V1
		DIVSS	XMM0, XMM1    					;questa ï¿½ la variabile rapporto rate/v
		SHUFPS	XMM0, XMM0, 00000000		
	
		MOVAPS	[rapportoB], XMM0

		
		MOV 	EDI, 0
		MOV 	ECX, 0
		MOV	EBX, [EBP + theta]
		
forj3B:	
		MOVAPS 	XMM2, [EAX + EDI]
		MOVAPS 	XMM3, [EAX + EDI + 16]
		MOVAPS 	XMM4, [EAX + EDI + 32]
		MOVAPS 	XMM5, [EAX + EDI + 48]
		
		MULPS 	XMM2, [rapportoB]
		MULPS 	XMM3, [rapportoB]
		MULPS 	XMM4, [rapportoB]
		MULPS 	XMM5, [rapportoB]
		
		MOVAPS 	XMM0, [EBX + EDI]
		MOVAPS 	XMM1, [EBX + EDI + 16]
		MOVAPS 	XMM6, [EBX + EDI + 32]
		MOVAPS 	XMM7, [EBX + EDI + 48]
		
		SUBPS 	XMM0, XMM2					;salvataggio in theta dei valori calcolati
		SUBPS 	XMM1, XMM3		
		SUBPS 	XMM6, XMM4		
		SUBPS 	XMM7, XMM5	
		
		MOVAPS	[EBX + EDI], XMM0
		MOVAPS 	[EBX + EDI + 16], XMM1
		MOVAPS 	[EBX + EDI + 32], XMM6
		MOVAPS 	[EBX + EDI + 48], XMM7
		
		;questo chiude il forj3
		ADD		ECX, 16
		ADD		EDI, elementiunrollingB
		CMP	ECX, [lengthB]
		JL		forj3B
		
		;ripristino gli indici v e i in modo da trovarli aggiornati al prossimo for
		MOV 	ECX, [vB]
		MOV 	EBX, [iB]
		
		; questo chiude il fori
		ADD 	EBX, [kB]
		CMP	EBX, [nB]
		JL 		foriB
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp								; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante
		
