%include "sseutils64.nasm"

section .data
	dim 		equ 8
	elementiunrolling 	equ	128
	align 32
	epsilon  	dq 0.00000001,0.00000001,0.00000001,0.00000001
	
section .bss
	;------------------------------------------->Variabili Adagrad
	alignb 32
	p		resq 1
	alignb 32
	i		resq 1
	alignb 32
	pfin		resq 1
	alignb 32
	lenght	resq 1	
	alignb 32
	theta	resq 1
	alignb 32
	tmp		resq 4
	alignb 32
	gj	  	resq 1
	alignb 32
	Gj	  	resq 1
	alignb 32
	sommatoria 	resq 1
	alignb 32
	osservazioni  	resq 1
	alignb 32
	y	  	resq 1
	alignb 32
	rapporto	resq 4
	alignb 32
	indice	resq 1
	alignb 32
	rate		resq 1
	alignb 32
	aux		resq 1
	alignb 32
	plenghtdim	resq 1
	
	;----------------------------------------------
	; 	VARIABILI AUSILIARIE
	;----------------------------------------------
	alignb 32
	auxxxxx 	resq 1
	alignb 32
	stampa	resq 4
	
	
section .text
	
global batch64Prod
global adagrad64Prod
global aggiornamentoG

batch64Prod:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	RBP				; salva il Base Pointer
		mov		RBP, RSP		; il Base Pointer punta al Record di Attivazione corrente
		pushaq					; salva i registri generali
		
		
		MOV	[p], RDI
		
		MOV	[pfin], RSI
		
		MOV	[i], RDX
		
		MOV	[osservazioni], RCX
		
		MOV	[lenght], R8
		
		MOV	[y], R9
		
		MOV	RAX, [RBP + 16]
		MOV	[tmp], RAX
		
		MOV	RAX, [RBP + 24]
		MOV	[theta], RAX
		
forpB:
		MOV	RDX, [p]
		
		VXORPD	YMM1, YMM1				;XMM1 contiene prodScal
		XOR 	RDI, RDI					;EDI indice j del for del prodotto scalare
		XOR 	RSI, RSI
		MOV 	RCX, RDX				;metto in ECX il valore di p
		IMUL	RCX, [lenght]				;moltiplico per lenght
		IMUL	RCX, RCX, dim				;IMUL ci salva p*lenght*dimB
		
		MOV	RBX,[theta]
		MOV	RDX,[osservazioni]	
		
forj1B:	
		VMOVAPD	YMM2, [RBX + RDI]
		VMULPD	YMM2, YMM2, [RDX + RCX]
		VADDPD	YMM1, YMM1, YMM2
		

		VMOVAPD	YMM3, [RBX + RDI + 32]
		VMULPD	YMM3, YMM3, [RDX + RCX + 32]
		VADDPD	YMM1, YMM1, YMM3
		
		
		VMOVAPD	YMM4, [RBX + RDI + 64]
		VMULPD	YMM4, YMM4, [RDX + RCX + 64]
		VADDPD	YMM1, YMM1, YMM4
		
		
		VMOVAPD	YMM5, [RBX + RDI + 96]
		VMULPD	YMM5, YMM5, [RDX + RCX + 96]
		VADDPD	YMM1, YMM1, YMM5	
		
		;questo chiude il forj1
		ADD		RSI, 16 
		ADD		RDI, elementiunrolling
		ADD 	RCX, elementiunrolling
		CMP	RSI, [lenght] 
		JL		forj1B

		VHADDPD	YMM1, YMM1, YMM1
		VHADDPD YMM1, YMM1, YMM1				;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		MOV	RBX, [y]
		MOV	RDX, [p]
		
		XOR		RDI, RDI
		IMUL	RDI, RDX, dim
		
		VSUBSD	XMM1, XMM1, [RBX + RDI]		;in XMM1 ho prodScal = prodScal - y[p]
		VMOVSD	[auxxxxx], XMM1
		VBROADCASTSD	YMM1, [auxxxxx]
			
		
		XOR		RDI,RDI					;EDI indice j del for del prodotto scalare
		MOV 	RSI, RDX
		SUB		RSI, [i]
		IMUL	RSI, [lenght]
		
		IMUL 	RSI, RSI, dim				;calcolo indice = p-i * lenght * dim
		
		
		
		IMUL	RDX, [lenght]				;posso riutilizzare EDX tanto il valore dell'indice p e salvato nella variabile p
		IMUL	RDX, RDX, dim			;IMUL ci salva p*lenght*dim in EDX
		
		MOV 	RAX, [tmp]
		MOV	RBX, [osservazioni]
		
		
forj2B:	
		VMOVAPD	YMM0, YMM1				
		VMOVAPD	YMM2, YMM1
		VMOVAPD	YMM4, YMM1
		VMOVAPD	YMM6, YMM1
			
		VMULPD	YMM0, YMM0, [RBX + RDX]
		VMULPD	YMM2, YMM2, [RBX + RDX + 32]
		VMULPD	YMM4, YMM4, [RBX + RDX + 64]
		VMULPD	YMM6, YMM6, [RBX + RDX + 96]		
		
		VMOVAPD		[RAX + RSI], XMM0		;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		VEXTRACTF128 XMM0,YMM0,1
		VMOVAPD		[RAX + RSI + 16], XMM0
		
		VMOVAPD		[RAX + RSI + 32], XMM2
		VEXTRACTF128 XMM2,YMM2,1
		VMOVAPD		[RAX + RSI + 48], XMM2
		
		VMOVAPD		[RAX + RSI + 64], XMM4
		VEXTRACTF128 XMM4,YMM4,1
		VMOVAPD		[RAX + RSI + 80], XMM4
		
		VMOVAPD		[RAX + RSI + 96], XMM6
		VEXTRACTF128 XMM6,YMM6,1
		VMOVAPD		[RAX + RSI + 112], XMM6
		
		;questo chiude il forj2
		ADD		RDI, 16
		ADD		RSI, elementiunrolling
		ADD		RDX, elementiunrolling
		CMP	RDI, [lenght]
		JL		forj2B
		
		; questo chiude il forp
		MOV	RDX, [p]
		INC 		RDX
		MOV	[p], RDX
		MOV 	RSI, [pfin]
		CMP 	RDX, RSI
		JL		forpB

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq
		mov	rsp, rbp								; ripristina lo Stack Pointer
		pop	rbp									; ripristina il Base Pointer
		ret
		
		
		
		
		
		
		
		
adagrad64Prod:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push	RBP				; salva il Base Pointer
		mov		RBP, RSP			; il Base Pointer punta al Record di Attivazione corrente
		pushaq		
		
		MOV	[p], RDI
		
		MOV	[pfin], RSI
		
		MOV	[theta], RDX
		
		MOV	[osservazioni], RCX
		
		MOV	[lenght], R8
		
		MOV	[gj], R9
		
		MOV	RAX, [RBP + 16]
		MOV	[Gj], RAX
		
		MOV	RAX, [RBP + 24]
		MOV	[y], RAX
		
		MOV	RAX, [RBP + 32]
		MOV	[i], RAX
		
		
		XOR		RAX, RAX
		XOR		RBX, RBX
		XOR		RCX, RCX
		XOR		RDX, RDX
		XOR		RDI, RDI
		XOR		RSI, RSI
		XOR		R8, R8
		XOR		R9, R9
		
		
forp:
		MOV 	RDX, [p]
		
		VXORPD	YMM1, YMM1				;XMM1 contiene prodScal
		XOR		RSI, RSI
		XOR 	RDI, RDI					;EDI indice j del for del prodotto scalare
		MOV 	ECX, EDX				;metto in ECX il valore di p
		IMUL	ECX, [lenght]				;moltiplico per lenght
		IMUL	ECX, ECX, dim				;IMUL ci salva p*lenght*dim
		
		MOV	[plenghtdim], ECX
		
		XOR		RBX, RBX
		XOR		RDX, RDX
		
		MOV	RBX, [theta]
		MOV	RDX, [osservazioni]
		
		ADD		RDX, [plenghtdim]
		
forj1:
		VMOVAPD	YMM2, [RBX]
		VMULPD	YMM2, YMM2, [RDX]
		VADDPD	YMM1, YMM1, YMM2
		VMOVAPD	YMM3, [RBX + 32]
		VMULPD	YMM3, YMM3, [RDX + 32]
		VADDPD	YMM1, YMM1, YMM3
		VMOVAPD	YMM4, [RBX + 64]
		VMULPD	YMM4, YMM4, [RDX + 64]
		VADDPD	YMM1, YMM1, YMM4
		VMOVAPD YMM5, [RBX + 96]
		VMULPD	YMM5, YMM5, [RDX + 96]
		VADDPD	YMM1, YMM1, YMM5	
		
		;questo chiude il forj1
		ADD		ESI, 16 
		ADD		RDX, elementiunrolling
		ADD 	RBX, elementiunrolling
		CMP	ESI, [lenght] 
		JL		forj1
		
	
		
		VHADDPD	YMM1, YMM1
		VHADDPD YMM1, YMM1			;nella prima posizione di XMM1 ho l'effettivo valore di prodScal
		
		XOR		R9, R9
		XOR 	RDX, RDX
		XOR		RBX, RBX
		
		
		MOV	R9, [y]
		
		MOV 	EDX, [p]
		IMUL	EDI, EDX, dim
		
		VSUBSD	XMM1, [R9 + RDI]		;in XMM1 ho prodScal = prodScal - y[p] 
		
		VMOVSD	[auxxxxx], XMM1
		
		VBROADCASTSD	YMM1, [auxxxxx]
		
		
		XOR		RDI, RDI					;EDI indice j del for del prodotto scalare
		XOR		RSI, RSI
		MOV 	ESI, EDX
		SUB		ESI, [i]
		IMUL	ESI, [lenght]
		IMUL 	ESI, ESI, dim				;calcolo indice = p-i * lenght * dim
		
		
		IMUL	EDX, [lenght]				;posso riutilizzare EDX tanto il valore dell'indice p e salvato nella variabile p
		IMUL	EDX, EDX, dim			;IMUL ci salva p*lenght*dim in EDX

		
		MOV 	RAX, [gj]
		MOV	RBX, [osservazioni]
		MOV	RCX, [Gj]
		
		ADD		RBX, RDX
		ADD		RAX, RSI
		ADD		RCX, RSI
		
		
forj2:	
		VMOVAPD	YMM0, YMM1				
		VMOVAPD	YMM2, YMM1
		VMOVAPD	YMM4, YMM1
		VMOVAPD	YMM6, YMM1
			
		VMULPD	YMM0, YMM0, [RBX]
		VMULPD	YMM2, YMM2, [RBX + 32]
		VMULPD	YMM4, YMM4, [RBX + 64]
		VMULPD	YMM6, YMM6, [RBX + 96]		
		
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
		
		
		VMOVAPD		[RCX], XMM0		;salvataggio in tmp dei valori calcolati, ovvero tutto ciò che sta a destra della sommatoria
		VEXTRACTF128 XMM0,YMM0,1
		VMOVAPD		[RCX + 16], XMM0
		
		VMOVAPD		[RCX + 32], XMM2
		VEXTRACTF128 XMM2,YMM2,1
		VMOVAPD		[RCX + 48], XMM2
		
		VMOVAPD		[RCX + 64], XMM4
		VEXTRACTF128 XMM4,YMM4,1
		VMOVAPD		[RCX + 80], XMM4
		
		VMOVAPD		[RCX + 96], XMM6
		VEXTRACTF128 XMM6,YMM6,1
		VMOVAPD		[RCX + 112], XMM6
		
		;questo chiude il forj2
		ADD		EDI, 16
		ADD		RAX, elementiunrolling
		ADD		RBX, elementiunrolling
		ADD		RCX, elementiunrolling
		CMP	EDI, [lenght]
		JL		forj2
		
		
		; questo chiude il forp
		XOR		RDX, RDX
		
		MOV	EDX, [p]
		INC 		EDX
		MOV	[p], EDX
		
		XOR		RSI, RSI
		MOV 	ESI, [pfin]
		CMP 	EDX, ESI
		JL		forp
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq
		mov	rsp, rbp								; ripristina lo Stack Pointer
		pop	rbp									; ripristina il Base Pointer
		ret