%include "sseutils32.nasm"

section .data
	dim 		equ 4
	unroll 	equ 4
	elementiunrolling 	equ	64
	align 16
	epsilon  	dd 0.00000001,0.00000001,0.00000001,0.00000001
	
section .bss
	;------------------------------------------->Variabili Adagrad
	alignb 16
	p		resd 1
	alignb 16
	i		resd 1
	alignb 16
	pfin		resd 1
	alignb 16
	lenght	resd 1	
	alignb 16
	theta	resd 1
	alignb 16
	tmp	resd 1
	alignb 16
	gj	  	resd 1
	alignb 16
	Gj	  	resd 1
	alignb 16
	sommatoria 	resd 1
	alignb 16
	osservazioni  	resd 1
	alignb 16
	y	  	resd 1
	alignb 16
	rapporto	resd 4
	alignb 16
	indice	resd 1
	alignb 16
	rate		resd 1
	alignb 16
	aux		resd 1
	
	;----------------------------------------------
	; 	VARIABILI AUSILIARIE
	;----------------------------------------------
	
	
section .text
	;---------------------------------------------indici ciclo prodScal Adagrad
	indicep			equ 8
	indicepfin			equ 12
	indicetheta 		equ 16
	indiceosservazioni	equ 20
	indicelunghezza	equ 24
	indicegj			equ 28
	indiceGj			equ 32
	indicey			equ 36
	indicei			equ 40
	;---------------------------------------------indici ciclo aggiornament sommatoria Adagrad
	Sindicep			equ 8
	Sindicepfin		equ 12
	Sindicelunghezza	equ 16
	Sindicesommatoria	equ 20
	Sindicerate		equ 24
	SindiceGj			equ 28
	Sindicegj			equ 32
	;---------------------------------------------indici ciclo prodScal batch
	Bindicep			equ 8
	Bindicepfin		equ 12
	Bindicetheta		equ 16
	Bindiceosservazioni 	equ 20
	Bindicelenght		equ 24
	Bindicetmp		equ 28
	Bindicey			equ 32
	Bindicei			equ 36
	
global adagrad32Prod
global aggiornamentoG
global batch32Prod

batch32Prod:
		; ------------------------------------------------------------
		;	 Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	ebp							; salva il Base Pointer
		mov		ebp, esp						; il Base Pointer punta al Record di Attivazione corrente
		push	ebx							; salva i registri da preservare
		push	esi
		push	edi	
		
		MOV 	EAX, [EBP + Bindicep]
		MOV	[p], EAX
		
		MOV 	EAX, [EBP + Bindicei]
		MOV	[i], EAX
		
		MOV 	EAX, [EBP + Bindicepfin]
		MOV	[pfin], EAX
		
		MOV	EAX, [EBP + Bindicetheta]
		MOV	[theta], EAX
		
		MOV	EAX, [EBP + Bindiceosservazioni]
		MOV	[osservazioni], EAX
		
		MOV	EAX, [EBP + Bindicelenght]
		MOV	[lenght], EAX
		
		MOV	EAX, [EBP + Bindicetmp]
		MOV	[tmp], EAX
		
		MOV	EAX, [EBP + Bindicey]
		MOV	[y], EAX
		
forpB:
		MOV	EDX, [p]
		
		XORPS	XMM1, XMM1				;XMM1 contiene prodScal
		XOR 	EDI, EDI					;EDI indice j del for del prodotto scalare
		XOR 	ESI, ESI
		MOV 	ECX, EDX				;metto in ECX il valore di p
		IMUL	ECX, [lenght]				;moltiplico per lenght
		IMUL	ECX, ECX, dim				;IMUL ci salva p*lenght*dimB
		
		MOV	EBX,[theta]
		MOV	EDX,[osservazioni]	
		
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
		ADD		EDI, elementiunrolling
		ADD 	ECX, elementiunrolling
		CMP	ESI, [lenght] 
		JL		forj1B

		HADDPS	XMM1, XMM1
		HADDPS	XMM1, XMM1				;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		MOV	EBX, [y]
		MOV	EDX, [p]
		
		XOR		EDI, EDI
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
		
		MOV 	EAX, [tmp]
		MOV	EBX, [osservazioni]
		
forj2B:	
		MOVAPS	XMM0, XMM1				
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
		
		;questo chiude il forj2
		ADD		EDI, 16
		ADD		ESI, elementiunrolling
		ADD		EDX, elementiunrolling
		CMP	EDI, [lenght]
		JL		forj2B
		
		; questo chiude il forp
		MOV	EDX, [p]
		INC 		EDX
		MOV	[p], EDX
		MOV 	ESI, [pfin]
		CMP 	EDX, ESI
		JL		forpB
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp								; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret	












adagrad32Prod:
		; ------------------------------------------------------------
		;	 Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	ebp							; salva il Base Pointer
		mov		ebp, esp						; il Base Pointer punta al Record di Attivazione corrente
		push	ebx							; salva i registri da preservare
		push	esi
		push	edi		
		
		
		MOV 	EAX, [EBP + indicep]
		MOV	[p], EAX
		
		MOV 	EAX, [EBP + indicei]
		MOV	[i], EAX
		
		MOV 	EAX, [EBP + indicepfin]
		MOV	[pfin], EAX
		
		MOV	EAX, [EBP + indicetheta]
		MOV	[theta], EAX
		
		MOV	EAX, [EBP + indiceosservazioni]
		MOV	[osservazioni], EAX
		
		MOV	EAX, [EBP + indicelunghezza]
		MOV	[lenght], EAX
		
		MOV	EAX, [EBP + indicegj]
		MOV	[gj], EAX
		
		MOV	EAX, [EBP + indiceGj]
		MOV	[Gj], EAX
		
		MOV	EAX, [EBP + indicey]
		MOV	[y], EAX
		
		
forp:
		MOV 	EDX, [p]
		
		XORPS	XMM1, XMM1				;XMM1 contiene prodScal
		XOR		ESI, ESI
		XOR 	EDI, EDI					;EDI indice j del for del prodotto scalare
		MOV 	ECX, EDX
		IMUL	ECX, [lenght]
		IMUL	ECX, ECX, dim
		
		MOV	EBX, [theta]
		MOV	EDX, [osservazioni]
		
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
		
		MOV	EBX, [y]
		MOV	EDX, [p]
		
		XOR		EDI, EDI
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
		
		MOV 	EAX, [gj]
		MOV	EBX, [osservazioni]
		MOV	ECX, [Gj]
		
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
		INC 		EDX
		MOV	[p], EDX
		MOV 	ESI, [pfin]
		CMP 	EDX, ESI
		JL		forp
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp								; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret	
		
	
	
	
	
	
	
	
	
	
aggiornamentoG:
		; ------------------------------------------------------------
		;	 Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	ebp							; salva il Base Pointer
		mov		ebp, esp						; il Base Pointer punta al Record di Attivazione corrente
		push	ebx							; salva i registri da preservare
		push	esi
		push	edi		
		
		MOV 	EAX, [EBP + Sindicep]
		MOV	[p], EAX
		MOV	[i], EAX
		
		MOV 	EAX, [EBP + Sindicepfin]
		MOV	[pfin], EAX
		
		MOV	EAX, [EBP + Sindicelunghezza]
		MOV	[lenght], EAX
		
		MOV	EAX, [EBP + Sindicesommatoria]
		MOV	[sommatoria], EAX
		
		MOV	EAX, [EBP + Sindicerate]
		MOV	[rate], EAX
		
		MOV	EAX, [EBP + SindiceGj]
		MOV	[Gj], EAX
		
		MOV	EAX, [EBP + Sindicegj]
		MOV	[gj], EAX
		
		MOV	ESI, [i]					;prendo i iniziale
		MOV	EDI, [pfin]				;metto la fine		
		
forp2:
		MOV 	EAX, [p]
		SUB		EAX, [i]
		IMUL	EAX, [lenght]
		IMUL 	EAX, EAX, dim				;calcolo indice = p-i * lenght * dim e lo metto in EAX
		MOV	[indice], EAX
		
		XOR		EBX, EBX
		MOV	ECX, [gj]
		
forj3:

		MOV	[indice], EAX
		MOV	EDX, [Gj]
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
		
		MOV	EDX, [sommatoria]
		
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
		INC		ESI
		MOV	[p],ESI	
		CMP	ESI, EDI
		JL		forp2
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp								; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret	















	
		
		
		
		
		
		
		
		
		
		
		
		
		