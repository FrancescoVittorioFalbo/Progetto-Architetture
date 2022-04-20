%include "sseutils64.nasm"

section .data
	dim 				equ 8
	elementiunrolling 	equ 128
	align 32
	epsilon  	dd 0.00000001,0.00000001,0.00000001,0.00000001

section .bss
	;-------------------------------------------------------------------------
	;		VARIABILI PER BATCH
	;-------------------------------------------------------------------------
	alignb 32
	kB		resd 1
	alignb 32
	nB		resd 1
	alignb 32
	lenghtB	resd 1	;numero di coordinate theta per ogni oggetto
	alignb 32
	osservB	resq 1	;indirizzo verso array osservazioni
	alignb 32
	yB		resq 1
	alignb 32
	thetaB	resq 1	;indirizzo verso array theta
	alignb 32
	indirizzotempB resq 1
	alignb 32
	tmpB	resq 4
	alignb 32
	tmp2B	resq 1
	alignb 32
	rateB 	resq 1
	alignb 32
	rapportoB	resq 4
	alignb 32
	plenghtdimB resq 1
	
	;----------------------------------------------
	; 	VARIABILI AUSILIARIE
	;----------------------------------------------
	
	alignb 32
	vB		resd 1
	alignb 32
	limiteB 	resd 1
	alignb 32
	iB		resd	1
	alignb 32		
	pB		resd	1
	alignb 32
	provaB	resq 1
	alignb 32
	prodScalB	resq 1
	alignb 32
	auxB		resq 1
	;-------------------------------------------------------------------------
	;		VARIABILI PER ADAGRAD
	;-------------------------------------------------------------------------

	alignb 32
	k		resq 1
	alignb 32
	n		resq 1
	alignb 32
	lenght	resq 1	;numero di coordinate theta per ogni oggetto
	alignb 32
	y		resq 1
	alignb 32
	rate 		resq	1
	alignb 32
	Gj		resq 1
	alignb 32
	gj		resq 1
	alignb 32
	theta	resq 1
	alignb 32
	osserv 	resq 1
	alignb 32
	sommatoria resq 1
	alignb 32
	plenghtdim	resq 1
	alignb 32
	prodScal	resq 1
	alignb 32
	rapporto	resq 4
	alignb 32
	indice	resq	1
	
	;----------------------------------------------
	; 	VARIABILI AUSILIARIE
	;----------------------------------------------
	
	alignb 32
	v		resq 4
	alignb 32
	limite 	resq 1	;contiene i+v
	alignb 32
	i		resq	1
	alignb 32		
	p		resq	1
	alignb 32
	aux		resq 1
	alignb 32
	prova	resq 1
	
	
section .text
	
; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro


global batch64
global adagrad64

batch64:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	RBP				; salva il Base Pointer
		mov		RBP, RSP		; il Base Pointer punta al Record di Attivazione corrente
		pushaq					; salva i registri generali				
		
		;PASSAGGIO DI LENGHT
		MOVSD		XMM0, [RCX]
		VEXTRACTPS 	EAX, XMM0,0
		MOV		[lenghtB], EAX
		
		
		;PASSAGGIO DI INPUT
		
		MOV	RAX, [RDI + 24]
		MOV	[nB], RAX
		
		MOV	RAX, [RDI + 32]
		MOV	[kB], RAX
		
		MOV   	RAX, [RDI+8]
		MOV	[yB], RAX
		
		VMOVSD   XMM7, [RDI+40]
		VMOVSD    [rateB], XMM7
		
		;PASSAGGIO DI THETA
		MOV	[thetaB], RDX
		
		MOV	[osservB], RSI
		
		XOR		RAX, RAX
		
		getmem	dim, [lenghtB]				;EAX non lo utilizziamo perche sara utilizzato da getmem
		
		MOV	[indirizzotempB], RAX
		
		XOR 	RBX, RBX					;EBX indice i del for dei passi di batch

foriB:
		
		MOV	[iB], EBX
		XOR		RSI, RSI					;ESI indirizzo di partenza dei 4 float di interesse
		XOR		RDI, RDI					;EDI numero di valori che scorriamo in un ciclo, (indice l del for di azzeramento di tmp)
		MOV	RAX, [indirizzotempB]
forlB:
		
		VXORPD 	YMM0, YMM0				;XMM0 ho tmp 
		VMOVAPD	[ RAX ], YMM0
		VMOVAPD	[ RAX + 32], YMM0
		VMOVAPD	[ RAX + 64], YMM0
		VMOVAPD	[ RAX + 96], YMM0
		ADD 	RAX, elementiunrolling
		ADD 	EDI, 16
		CMP 	EDI, [lenghtB]
		JL 		forlB
		
	
		XOR		RCX, RCX
		;----------------------
		;CALCOLO DI V
		;---------------------
		MOV 	ECX, [nB] 					;In ECX inserisco il valore v, contatore di elementi nel corrente passo di batch
		SUB 	ECX, EBX
		CMP 	ECX, [kB]
		JLE 		saltoB
		
		MOV 	ECX, [kB]
saltoB:
		VCVTSI2SD	XMM0, ECX				;convertire integer to float
		MOVSD	[vB], XMM0
		;MOV	[v], ECX
		
		MOV 	EDX, EBX 				;EDX indice p per scorrere gli elementi del corrente passo di batch faccio p=i
		MOV 	ESI, EBX
		ADD 	ESI, ECX 					;In ESI ho inserito i+v per il forp
										;edx ho p e in esi ho i+v condizione di inizio del forp
		
		MOV	[limiteB], ESI

forpB:
		VXORPD	YMM1, YMM1				;XMM1 contiene prodScal
		MOV 	EDI, 0					;EDI indice j del for del prodotto scalare
		MOV 	ESI, 0
		
		XOR		RCX, RCX
		
		MOV 	ECX, EDX				;metto in ECX il valore di p
		IMUL	ECX, [lenghtB]				;moltiplico per lenght
		IMUL	ECX, ECX, dim				;IMUL ci salva p*lenght*dim
		
		
		MOV 	[pB], EDX
		MOV	RBX,[thetaB]
		MOV	RDX,[osservB]
		
		
		MOV	[plenghtdimB], ECX
		ADD		RDX, [plenghtdimB]
forjB1:	
		VMOVAPD	YMM2, [RBX]
		VMULPD	YMM2, YMM2, [RDX]
		VADDPD	YMM1, YMM2
		VMOVAPD	YMM3, [RBX + 32]
		VMULPD	YMM3, YMM3, [RDX + 32]
		VADDPD	YMM1, YMM3
		VMOVAPD	YMM4, [RBX+ 64]
		VMULPD	YMM4, YMM4, [RDX + 64]
		VADDPD	YMM1, YMM4
		VMOVAPD	YMM5, [RBX+ 96]
		VMULPD	YMM5, YMM5, [RDX + 96]
		VADDPD	YMM1, YMM5	
		
		;questo chiude il forj1
		ADD		ESI, 16 
		ADD		RDX, elementiunrolling
		ADD		RBX, elementiunrolling
		CMP	ESI, [lenghtB] 
		JL		forjB1
		
		
		VHADDPD	YMM1, YMM1
		VHADDPD YMM1, YMM1				;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		XOR		R9, R9
		MOV	R9, [yB]
		
		XOR 	RDX, RDX
		XOR		RBX, RBX
		
		MOV 	EDX, [pB]
		
		IMUL	EDI, EDX, dim
		
		MOV	[auxB], EDI
		ADD		R9, [auxB]
		
		VSUBSD	XMM1, [R9]		;in XMM1 ho prodScal = prodScal - y[p] 
		VMOVSD	[prodScalB], XMM1
		
		MOV 	EDI, 0					;EDI indice j del for del prodotto scalare
		IMUL	EDX, [lenghtB]				;posso riutilizzare EDX tanto il valore dell'indice p e salvato nella variabile p
		IMUL	EDX, EDX, dim			;IMUL ci salva p*lenght*dim
		
		MOV	RBX, [osservB]
		ADD		RBX, RDX
	
	
		MOV	RAX, [indirizzotempB]
		
forjB2:	
		VBROADCASTSD	YMM0, [prodScalB]
		VMULPD			YMM0, [RBX]
		
		VBROADCASTSD	YMM2, [prodScalB]
		VMULPD			YMM2, [RBX + 32]	
		
		VBROADCASTSD	YMM4, [prodScalB]
		VMULPD			YMM4, [RBX + 64]
		
		VBROADCASTSD	YMM6, [prodScalB]
		VMULPD			YMM6, [RBX +  96]
		
		VADDPD	YMM0, [RAX]			;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		VADDPD	YMM2, [RAX + 32]
		VADDPD	YMM4, [RAX + 64]
		VADDPD	YMM6, [RAX + 96]
		
		VMOVAPD	[RAX], YMM0			;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		VMOVAPD	[RAX + 32], YMM2
		VMOVAPD	[RAX + 64], YMM4
		VMOVAPD	[RAX + 96], YMM6
		
		;questo chiude il forj2
		ADD		EDI, 16
		ADD		RAX, elementiunrolling
		ADD		RBX, elementiunrolling
		CMP	EDI, [lenghtB]
		JL		forjB2
		
		XOR		RBX, RBX
		
		; questo chiude il forp
		XOR		RDX, RDX
		MOV	EDX, [pB]
		MOV 	ESI, [limiteB]
		ADD		EDX, 1
		CMP 	EDX, ESI
		JL		forpB
		
		
		VMOVSD			XMM0, [rateB]	
		VMOVSD			XMM1,[vB]			;FACCIAMO DUE CONVERSIONI
		VDIVSD			XMM0, XMM1    					;questa è la variabile rapporto rate/v
	
		VMOVSD			[rapportoB], XMM0
		VBROADCASTSD	YMM0, [rapportoB]		
		VMOVAPD			[rapportoB], YMM0
		
		
		MOV 	ECX, 0
		
		
		MOV	RBX, [thetaB]
		MOV	RAX, [indirizzotempB]
forjB3:	
		VMOVAPD 	YMM2, [RAX]
		VMOVAPD 	YMM3, [RAX + 32]
		VMOVAPD 	YMM4, [RAX + 64]
		VMOVAPD 	YMM5, [RAX + 96]
		
		VMULPD 		YMM2, YMM2, [rapportoB]
		VMULPD		YMM3, YMM3, [rapportoB]
		VMULPD		YMM4, YMM4, [rapportoB]
		VMULPD		YMM5, YMM5, [rapportoB]
		
		VMOVAPD 	YMM0, [RBX]
		VMOVAPD 	YMM1, [RBX + 32]
		VMOVAPD 	YMM6, [RBX + 64]
		VMOVAPD 	YMM7, [RBX + 96]
		
		VSUBPD 		YMM0, YMM2					;salvataggio in theta dei valori calcolati
		VSUBPD		YMM1, YMM3		
		VSUBPD		YMM6, YMM4		
		VSUBPD		YMM7, YMM5	
		
		VMOVAPD		[RBX], YMM0
		VMOVAPD 	[RBX + 32], YMM1
		VMOVAPD 	[RBX + 64], YMM6
		VMOVAPD 	[RBX + 96], YMM7
		
		;questo chiude il forj3
		ADD		ECX, 16
		ADD		RAX, elementiunrolling
		ADD		RBX, elementiunrolling
		CMP	ECX, [lenghtB]
		JL		forjB3		

		
		;ripristino gli indici v e i in modo da trovarli aggiornati al prossimo for
		XOR		RBX, RBX
		MOV 	EBX, [iB]
		
		; questo chiude il fori
		ADD 	EBX, [kB]
		CMP	EBX, [nB]
		JL 		foriB
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante
		








adagrad64:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	RBP				; salva il Base Pointer
		mov		RBP, RSP		; il Base Pointer punta al Record di Attivazione corrente
		pushaq					; salva i registri generali	
		
		
		;PASSAGGIO DI LENGHT		
		VCVTSD2SI	EAX, XMM0
		MOV		[lenght], EAX
	
		MOV	RAX, [RDI + 24]
		MOV	[n], RAX
		
		MOV	RAX, [RDI + 32]
		MOV	[k], RAX
		
		MOV   	RAX, [RDI+8]
		MOV	[y], RAX
		
		VMOVSD   XMM7, [RDI+40]
		VMOVSD    [rate], XMM7	
		
		MOV	[osserv], RSI
		
		MOV	[theta], RDX
		
		MOV	[Gj], RCX
		
		MOV	[gj], R8
		
		MOV	[sommatoria], R9
		
		XOR		RAX, RAX
		XOR		RBX, RBX
		XOR		RCX, RCX
		XOR		RDX, RDX
		XOR		RDI, RDI
		XOR		RSI, RSI
		XOR		R8, R8
		XOR		R9, R9
		
		
fori:
		MOV	[i], EBX
		MOV 	ECX, [n]		 ;In ECX inserisco il valore v, contatore di elementi nel corrente passo di batch
		SUB 	ECX, EBX
		CMP 	ECX, [k]
		
		JLE		salto
		
		MOV 	ECX, [k]
salto:
		VCVTSI2SD			XMM0, ECX		;convertire integer to float
		VMOVSD				[aux], XMM0
		VBROADCASTSD		YMM1, [aux]
		VMOVAPD				[v], YMM1
		
		MOV 	EDX, EBX ;EDX indice p per scorrere gli elementi del corrente passo di batch faccio p=i
		MOV 	ESI, EBX
		ADD 	ESI, ECX 	;In ESI ho inserito i+v per il forp
						;edx ho p e in esi ho i+v condizione di inizio del forp
		
		MOV	[limite], ESI
		
forp:
		MOV 	[p], EDX
		
		VXORPD	YMM1, YMM1				;YMM1 contiene prodScal
		XOR		RDI, RDI
		XOR		RSI, RSI
		MOV 	ECX, EDX				;metto in ECX il valore di p
		IMUL	ECX, [lenght]				;moltiplico per lenght
		IMUL	ECX, ECX, dim				;IMUL ci salva p*lenght*dim
		
		MOV	[plenghtdim], ECX
		
		MOV	RBX,[theta]
		MOV	RDX,[osserv]
		
		ADD		RDX, [plenghtdim]
		
forj1:	
		VMOVAPD	YMM2, [RBX]
		VMULPD	YMM2, [RDX]
		VADDPD	YMM1, YMM2
		VMOVAPD	YMM3, [RBX + 32]
		VMULPD	YMM3, [RDX + 32]
		VADDPD	YMM1, YMM3
		VMOVAPD	YMM4, [RBX + 64]
		VMULPD	YMM4, [RDX + 64]
		VADDPD	YMM1, YMM4
		VMOVAPD	YMM5, [RBX + 96]
		VMULPD	YMM5, [RDX + 96]
		VADDPD	YMM1, YMM5	
		
		;questo chiude il forj1
		ADD		ESI, 16 
		ADD		RDX, elementiunrolling
		ADD		RBX, elementiunrolling
		CMP	ESI, [lenght] 
		JL		forj1
		

		VHADDPD	YMM1, YMM1
		VHADDPD	YMM1, YMM1				;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		
		XOR		R9, R9
		XOR 	RDX, RDX
		XOR		RBX, RBX
		
		MOV	R9, [y]
		
		MOV 	EDX, [p]
		IMUL	EDI, EDX, dim

		;-------------------------------------------
		;ho inserito RDI cosi posso fare 
		;la somma tra due registri essendo a 64 bit 
		;il valore altrimenti non potevo
		;-------------------------------------------
		VSUBSD	XMM1, [R9 + RDI]		;in XMM1 ho prodScal = prodScal - y[p] 
		
		VMOVSD	[prodScal], XMM1
		VBROADCASTSD	YMM1, [prodScal]
		
		XOR		RDI, RDI					;EDI indice j del for del prodotto scalare
		XOR		RSI, RSI
		MOV 	ESI, EDX
		SUB		ESI, [i]
		IMUL	ESI, [lenght]
		IMUL 	ESI, ESI, dim				;calcolo indice = p-i * lenght * dim
		
		
		IMUL	EDX, [lenght]				;posso riutilizzare EDX tanto il valore dell'indice p e salvato nella variabile p
		IMUL	EDX, EDX, dim			;IMUL ci salva p*lenght*dim in EDX

		
		MOV 	RAX, [gj]
		MOV	RBX, [osserv]
		MOV	RCX, [Gj]
		
		ADD		RBX, RDX
		ADD		RAX, RSI
		ADD		RCX, RSI
		
forj2:	
		VMOVAPD	YMM0, YMM1		;in aux ce prodScal
		VMOVAPD	YMM2, YMM1
		VMOVAPD	YMM4, YMM1
		VMOVAPD	YMM6, YMM1
			
		VMULPD 	YMM0, [RBX]
		VMULPD 	YMM2, [RBX + 32]
		VMULPD	YMM4, [RBX + 64]
		VMULPD	YMM6, [RBX + 96]		
		
		VMOVAPD	[RAX], YMM0		;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		VMOVAPD	[RAX + 32], YMM2
		VMOVAPD	[RAX + 64], YMM4
		VMOVAPD	[RAX + 96], YMM6
		
		VMULPD	YMM0, YMM0
		VMULPD	YMM2, YMM2
		VMULPD	YMM4, YMM4
		VMULPD	YMM6, YMM6
		
		VADDPD	YMM0, [RCX]
		VADDPD	YMM2, [RCX + 32]
		VADDPD	YMM4, [RCX + 64]
		VADDPD	YMM6, [RCX + 96]
		
		VMOVAPD	[RCX ], YMM0		;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		VMOVAPD	[RCX + 32], YMM2
		VMOVAPD	[RCX + 64], YMM4
		VMOVAPD	[RCX + 96], YMM6
		
		;questo chiude il forj2
		ADD		EDI, 16
		ADD		RAX, elementiunrolling
		ADD		RBX, elementiunrolling
		ADD		RCX, elementiunrolling
		CMP	EDI, [lenght]
		JL		forj2
		
		; questo chiude il forp
		XOR		RDX, RDX
		XOR		RSI, RSI
		MOV	EDX, [p]
		MOV 	ESI, [limite]
		ADD		EDX, 1
		CMP 	EDX, ESI
		JL		forp

		XOR		RSI, RSI					;ESI indirizzo di partenza dei 4 float di interesse
		XOR		RDI, RDI				;EDI numero di valori che scorriamo in un ciclo, (indice l del for di azzeramento di tmp)
		MOV	RAX, [sommatoria]
forl:
		
		VXORPD 	YMM0, YMM0				;XMM0 ho sommatoria 
		VMOVAPD	[ RAX ], YMM0
		VMOVAPD	[ RAX + 32], YMM0
		VMOVAPD	[ RAX + 64], YMM0
		VMOVAPD	[ RAX + 96], YMM0
		ADD 	RAX, elementiunrolling			
		ADD 	EDI, 16
		CMP 	EDI, [lenght]
		JL 		forl
		
		XOR 	RSI, RSI
		MOV	ESI, [i]
		MOV	[p], ESI
		
forp2:
		XOR		RAX, RAX
		MOV 	EAX, [p]
		SUB		EAX, [i]
		IMUL	EAX, [lenght]
		IMUL 	EAX, EAX, dim				;calcolo indice = p-i * lenght * dim e lo metto in EAX
		
		MOV	[indice], EAX
	
		MOV	RCX, [gj]
		MOV	RDX, [Gj]
		
		ADD		RDX, [indice]
		ADD		RCX, [indice]
		
		MOV	RAX, [sommatoria]
		
		XOR		RSI, RSI
forj3:


		VMOVAPD	YMM4, [RDX]			;qui abbiamo Gj
		VMOVAPD	YMM5, [RDX + 32]
		VMOVAPD	YMM6, [RDX + 64]
		VMOVAPD	YMM7, [RDX + 96]
		
		VADDPD	YMM4, [epsilon]
		VADDPD	YMM5, [epsilon]
		VADDPD	YMM6, [epsilon]
		VADDPD	YMM7, [epsilon]
	
		;JMP return
		
		VSQRTPD	YMM4, YMM4
		VSQRTPD	YMM5, YMM5
		VSQRTPD	YMM6, YMM6
		VSQRTPD	YMM7, YMM7
		
		VBROADCASTSD YMM0, [rate]
		
		VMOVAPD 	YMM1, YMM0
		VMOVAPD		YMM2, YMM0
		VMOVAPD		YMM3, YMM0
		
		VDIVPD	YMM0, YMM4
		VDIVPD	YMM1, YMM5
		VDIVPD	YMM2, YMM6
		VDIVPD	YMM3, YMM7
		
		VMOVAPD	YMM4, [RCX]			;qui abbiamo gj
		VMOVAPD	YMM5, [RCX + 32]
		VMOVAPD	YMM6, [RCX + 64]
		VMOVAPD	YMM7, [RCX + 96]
		
	
		VMULPD	YMM0, YMM4
		VMULPD	YMM1, YMM5
		VMULPD	YMM2, YMM6
		VMULPD	YMM3, YMM7
	
		
		VADDPD	YMM0, [RAX]			;qui abbiamo sommatoria
		VADDPD	YMM1, [RAX + 32]
		VADDPD	YMM2, [RAX + 64]
		VADDPD	YMM3, [RAX + 96]
		
		VMOVAPD	[RAX], YMM0
		VMOVAPD	[RAX + 32], YMM1
		VMOVAPD	[RAX + 64], YMM2
		VMOVAPD	[RAX + 96], YMM3
		
		;questo chiude forj3
		ADD		ESI, 16
		ADD		RAX, elementiunrolling
		ADD		RCX, elementiunrolling
		ADD		RDX, elementiunrolling
		CMP	ESI, [lenght]
		JL		forj3
		
		;questo chiude forp2
		XOR		RDI, RDI
		MOV	EDI, [limite]
		XOR		RSI, RSI
		MOV	ESI, [p]
		ADD		ESI, 1
		MOV	[p],ESI	
		CMP	ESI, EDI
		JL		forp2
		
		
		XOR		RDI, RDI
		XOR		RCX, RCX
		XOR		RAX, RAX
		XOR		RBX, RBX
		
		
		MOV	RBX, [theta]
		MOV	RAX, [sommatoria]
		
forj4:	
		VMOVAPD YMM2, [RAX]
		VMOVAPD YMM3, [RAX + 32]
		VMOVAPD YMM4, [RAX + 64]
		VMOVAPD YMM5, [RAX + 96]
		
		VDIVPD	YMM2, [v]
		VDIVPD	YMM3, [v]
		VDIVPD	YMM4, [v]
		VDIVPD	YMM5, [v]
		
		VMOVAPD  YMM0, [RBX]
		VMOVAPD  YMM1, [RBX + 32]
		VMOVAPD  YMM6, [RBX + 64]
		VMOVAPD  YMM7, [RBX + 96]
		
		VSUBPD	YMM0, YMM2					;salvataggio in theta dei valori calcolati
		VSUBPD	YMM1, YMM3		
		VSUBPD	YMM6, YMM4		
		VSUBPD	YMM7, YMM5	
		
		VMOVAPD  [RBX], YMM0
		VMOVAPD  [RBX + 32], YMM1
		VMOVAPD  [RBX + 64], YMM6
		VMOVAPD  [RBX + 96], YMM7
		
		;questo chiude il forj4
		ADD		ECX, 16
		ADD		RAX, elementiunrolling
		ADD		RBX, elementiunrolling
		CMP	ECX, [lenght]
		JL		forj4
		
		;ripristino gli indici v e i in modo da trovarli aggiornati al prossimo for
		XOR		RBX, RBX
		MOV 	EBX, [i]
		
		; questo chiude il fori
		ADD 	EBX, [k]
		CMP	EBX, [n]
		JL 		fori

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq						; ripristina i registri generali
		mov		rsp, rbp				; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante