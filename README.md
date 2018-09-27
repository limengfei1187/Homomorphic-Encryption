# Homomorphic-Encryption
This repository is about the implementation of homomorphic encryption scheme, include DGHV, CMNT and CNT. As the some time, we test efficiency of different homomorphic encryption scheme. We also show efficient of homomorphic evaluation of the integer arithmetic operations base on DGHV, CMNT and CNT schemes. 

Our code is built on top of GnuMP. We tested our implementation on a desktop computer with Intel Core i5-3470 running at 3.2GHz, on which we run an Ubuntu 18.04 with 8GB of RAM and with the gcc compiler version 6.2.

If you want use our program, you should ensure that GUN MP Library has been installed in your system. You can visit the GNU MP web pages at https://gmplib.org/ to download GUN MP Library.

In /Data/Data1.cvs, we show CPU-time of homomorphic evaluation of the integer arithmetic operations, include 
  Homomorphic Evaluation of the Complement Operations(HE-com),
  Homomorphic Evaluation of the Addition Operations(HE-add),
  Homomorphic Evaluation of the Subtraction Operations(HE-sub),
  Homomorphic Evaluation of the Multiplication Operations(HE-mul),
  Homomorphic Evaluation of the Division Operation(HE-div),
based on DGHV, CMNT and CNT schemes. We use the ciphertext vector of 16 and 8 bit plaintext to test the efficient of homomorphic evaluation of the integer arithmetic operations in different parameters. We show parameters include  L0, L1, L2, L3 and L4. As shown /Data/parameters.cvs, where λ, ρ^', η, γ and τ are importent parameters in  homomorphic encryption scheme.

In /Data/Data2.cvs, we show CPU-time of DGHV and CMNT scheme, include CPU-time of executing homomorphic addition (for single ciphertext), CPU-time of executing homomorphic multiplication (for single ciphertext), CPU-time of executing ciphertext-expending and ciphertext-refreshing in diferent parameters. Also, that parameters include  L0, L1, L2, L3 and L4.

In /DGHV/ directory, the file of "DGHV.c" is source code about DGHV scheme. You can use follow command to compile it in Linux.
  gcc DGHV.c -lgmp -lm -o DGHV
The file of /DGHV/DGHV was already compiled. You can run it in terminal of linux.
The file of /DGHV/privatekey storage private key.
The public key files is too large to upload it.

In /CMNT/ directory, the file of "DGHV.c" is source code about DGHV scheme. You can use follow command to compile it in Linux.
  gcc CMNT.c -lgmp -lm -o CMNT
The file of /CMNT/CMNT was already compiled. You can run it in terminal of linux.
The file of /CMNT/privatekey storage private key.
The public key files is too large to upload it.

In /CNT/ directory, the file of "DGHV.c" is source code about DGHV scheme. You can use follow command to compile it in Linux.
  gcc CNT.c -lgmp -lm -o CNT
The file of /CNT/CNT was already compiled. You can run it in terminal of linux.
The file of /CNT/privatekey storage private key.
The file of /CNT/publickey storage private key.


