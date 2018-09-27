/*************************************************************************
    > File Name: DGHV.c
    > Author: Mengfei Li
    > Mail: 1187337437@qq.com 
    > Created Time: 2018-09-27 Thursday 11:52:23
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>

#define L0   0
#define L1   1
#define L2   2
#define L3   3
#define L4   4


typedef struct paramters{
	unsigned long lam;     //安全参数
	unsigned long Rho;     //加密噪音
	unsigned long rho;     //公钥噪音
	unsigned long eat;     //私钥长度
	unsigned long gam;     //公钥长度
	unsigned long Theta;   //集合维度
	unsigned long theta;   //稀疏子集长度
	unsigned long tau;     //公钥个数
	unsigned long n;       //精度
}paramters;

double duration;
clock_t start,finish;
paramters paraSet;
int s[10000];
mpz_t pk[10000];
mpz_t p;
mpf_t y[10000];

mpz_t evla_matrix[10000][20];
mpz_t hw_matrix[10][10000];

void DGHV_init(const int para){
	
	if(para==L0){

		paraSet.lam = 52;
		paraSet.Rho = 24;
		paraSet.rho = 24;
		paraSet.eat = 1632;
		paraSet.gam = 2000000;	
		paraSet.Theta = 500;
		paraSet.theta = 15;
		paraSet.n = 4;
		paraSet.tau = 1000;
	}
	if(para==L1){

		paraSet.lam = 52;
		paraSet.Rho = 24;
		paraSet.rho = 24;
		paraSet.eat = 1728;
		paraSet.gam = 2500000;	
		paraSet.Theta = 500;
		paraSet.theta = 15;
		paraSet.n = 4;
		paraSet.tau = 1000;
	}
	if(para==L2){

		paraSet.lam = 52;
		paraSet.Rho = 24;
		paraSet.rho = 24;
		paraSet.eat = 1801;
		paraSet.gam = 3000000;	
		paraSet.Theta = 500;
		paraSet.theta = 15;
		paraSet.n = 4;
		paraSet.tau = 1000;
	}
	if(para==L3){

		paraSet.lam = 52;
		paraSet.Rho = 24;
		paraSet.rho = 24;
		paraSet.eat = 1874;
		paraSet.gam = 3500000;	
		paraSet.Theta = 500;
		paraSet.theta = 15;
		paraSet.n = 4;
		paraSet.tau = 1000;
	}
	if(para==L4){

		paraSet.lam = 52;
		paraSet.Rho = 24;
		paraSet.rho = 24;
		paraSet.eat = 1993;
		paraSet.gam = 4200000;	
		paraSet.Theta = 500;
		paraSet.theta = 15;
		paraSet.n = 4;
		paraSet.tau = 1000;
	}

	int i,j;
	for(i=0;i<=paraSet.Theta;i++){mpz_init_set_ui(hw_matrix[0][i],1);}
	for(i=0;i<8;i++){mpz_init_set_ui(hw_matrix[i][0],0);}
	
	mpf_set_default_prec(paraSet.gam+64);
	printf("%ld\n",mpf_get_default_prec () );
}


void gen_rrandomb(mpz_t rand_num,mp_bitcnt_t n){
	unsigned long seed;
	gmp_randstate_t rstate;

	struct timeval t_val;
	gettimeofday(&t_val, NULL);
	seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;

	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);

	mpz_rrandomb(rand_num,rstate,n);

}

void gen_urandomm(mpz_t rand_num,mpz_t upper_b){
	unsigned long seed;
	gmp_randstate_t rstate;

	struct timeval t_val;
	gettimeofday(&t_val, NULL);

	seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);

	mpz_urandomm(rand_num,rstate,upper_b);

}

void gen_prime_num(mpz_t prime_num,mp_bitcnt_t n){

	mpz_t rand_num;
	mpz_init(rand_num);
	gen_rrandomb(rand_num,n);
	mpz_nextprime(prime_num,rand_num);
	mpz_clear(rand_num);
	
}

void DGHV_encrypt(mpz_t ciphertext, unsigned long plaintext, unsigned long tau,mp_bitcnt_t Rho){
	
	int i,j,r;
	
	unsigned long seed;
	struct timeval t_val;
	mpz_t rand_num;
	
	mpz_init(rand_num);
	gettimeofday(&t_val, NULL);
	seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;
	srand(seed);
	for(i=0;i<88;i++){
		r=rand()%tau;
		while(r==0){r=rand()%tau;}
		mpz_add(ciphertext,ciphertext,pk[r]);
	}
	mpz_mul_ui(ciphertext,ciphertext,2);
	mpz_mod(ciphertext,ciphertext,pk[0]);
	gen_rrandomb(rand_num,Rho);
	mpz_mul_ui(rand_num,rand_num,2);
	mpz_add_ui(rand_num,rand_num,plaintext);
	mpz_add(ciphertext,ciphertext,rand_num);
	
	mpz_clear(rand_num);
}

unsigned long int DGHV_decrypt(mpz_t ciphertext){
	mpz_t plaintext;
	
	mpz_init(plaintext);
	mpz_mod(plaintext,ciphertext,p);
	mpz_mod_ui(plaintext,plaintext,2);
	return mpz_get_ui(plaintext);
	
	
}

void DGHV_add(mpz_t sum,mpz_t c1,mpz_t c2){
	mpz_add(sum,c1,c2);
	mpz_mod(sum,sum,pk[0]);
}

void DGHV_mul(mpz_t product,mpz_t c1,mpz_t c2){
	mpz_mul(product,c1,c2);
	mpz_mod(product,product,pk[0]);
}

void gen_pk_set(mp_bitcnt_t rho, mp_bitcnt_t gam, unsigned long tau){

	int i;
	mpz_t q,r,base,upper_b;
	mpz_init(q);
	mpz_init(r);
	mpz_init(upper_b);
	mpz_init_set_ui(base,2);

	mpz_pow_ui(upper_b,base,gam);
	mpz_fdiv_q(upper_b,upper_b,p);

	for(i=0;i<tau;i++){

		mpz_init(pk[i]);
		gen_urandomm(q,upper_b);
		gen_rrandomb(r,rho);
		if(i==0){
			if(mpz_odd_p(q)==0){
				mpz_sub_ui(q,q,1);
			}
			mpz_fdiv_q_ui(upper_b,q,4);
			mpz_mul(pk[i],p,q);
		}else{
			mpz_mul(pk[i],p,q);
			mpz_add(pk[i],pk[i],r);
		}
	}
	 
	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(upper_b);
	mpz_clear(base);
}

void div_trun(mpz_t q, mpz_t  n, mpz_t d){

	mpz_t r,s;

	mpz_init(r);
	mpz_init(s);

	mpz_fdiv_qr(q,r,n,d);
	mpz_fdiv_q_ui(s,d,2);

	if(mpz_cmp(r,s)>=0){
		mpz_add_ui(q,q,1);
	}

	mpz_clear(s);
	mpz_clear(r);
}

void expend_p2y(mp_bitcnt_t gam,unsigned long theta,unsigned long Theta){

	int i,j;
	unsigned long int seed,r;
	
	
	mpz_t yy[15];
	mpz_t quotient, remainder,result;
	mpz_t xp;
	mpz_t rand_num,ui;
	mpf_t numerator,denominator,base;

	mpz_init(quotient);
	mpz_init(remainder);
	mpz_init(result);
	mpz_init(xp);
	mpz_init(rand_num);
	mpz_init(ui);
	
	mpf_init(numerator);
	mpf_init(denominator);
	mpf_init_set_ui(base,2);

	for(i=0;i<Theta;i++){ 
		s[i]=0;
	}
	
	seed = time(NULL);
	srand(seed);
        for(i=0;i<theta;i++){
	    r= rand()%Theta;
	    //printf("%ld\n",r);
	    if(s[r]==0){
                s[r]=1;
		mpz_init_set_ui(yy[i],0);
	    }else if(s[r]==1){
		i=i-1;
	    }
	    //printf("%d\n",i);
		
	}

	mpz_ui_pow_ui(xp,2,gam+2);
	div_trun(xp,xp,p);
	mpz_fdiv_qr_ui(quotient,remainder,xp,theta);

	for(i=0; i<theta; i++){
		mpz_add(yy[i], yy[i], quotient);
	}
	mpz_add(yy[theta-1], yy[theta-1], remainder);
	
	for(i=0; i<theta; i++){
	
		mpz_fdiv_qr_ui(quotient,remainder,yy[i],rand()%Theta+1);
		mpz_add(result,quotient,remainder);
		mpz_sub(yy[i],yy[i],result);
		
		mpz_fdiv_qr_ui(quotient,remainder,result,theta);
		mpz_add(result,remainder,quotient);
		
		for(j=0;j<theta;j++){
			if(j==i){
				mpz_add(yy[j],yy[j],result);
				
			}else{
				mpz_add(yy[j],yy[j],quotient);
			}
		}
	}
	
	mpz_ui_pow_ui(ui,2,gam+3);
	mpf_pow_ui(denominator,base,gam+2);

	for(i=0,j=0;i<Theta;i++){
		mpf_init(y[i]);
		if(s[i]==0){
			gen_urandomm(rand_num,ui);
			mpf_set_z(numerator,rand_num);

		}else if(s[i]==1){
			mpf_set_z(numerator,yy[j]);
			j++;
		}
		mpf_div(y[i],numerator,denominator);
		//printf("%d\n",i);
	}

	mpz_clear(rand_num);
	mpz_clear(ui);
	mpz_clear(quotient);
	mpz_clear(remainder);
	mpz_clear(result);
	mpz_clear(xp);
	mpf_clear(base);
	mpf_clear(numerator);
	mpf_clear(denominator);
}

void save(const char* type){
    int i;
    FILE *out;
    if(strcmp(type,"publickey")==0){
        if((out = fopen("publickey","wt"))== NULL){
            fprintf(stderr,"Cannot open publickey file\n");
        }
        gmp_fprintf(out,"**********|pk_vector|**********\n");
        for(i=0;i<paraSet.tau;i++){
            gmp_fprintf(out,"%Zd\n",pk[i]);
        }
       
    }
    if(strcmp(type,"privatekey")==0){
        if((out = fopen("privatekey","wt"))== NULL){
            fprintf(stderr,"Cannot open publickey file\n");
        }
        gmp_fprintf(out,"**********|privatekey.p|**********\n");
        gmp_fprintf(out,"%Zd\n",p);
    }
   
    fclose(out);

}

void DGHV_genkey(paramters para){

	int i;
	double total=0;
	printf("The time of DGHV_genkey(paraSet) function start tunning...\n");

	printf("Generate sk:p\n");
	start = clock();
	gen_prime_num(p,para.eat);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
	printf("Generate success! The time:%f seconds\n",duration);
	
	printf("privateKey:\n");
	gmp_printf("%Zd\n",p);
	
	printf("Save private key...\n");
	save("privatekey");
	printf("Save private key end\n");
	
	printf("Generate pk set...\n");
	start = clock();
	gen_pk_set(para.rho,para.gam,para.tau);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
	printf("Generate pk set success,total:%ld, the time: %f seconds\n",para.tau,duration);
	
	printf("Save publick key...\n");
	save("publickey");
	printf("Save publick key end\n");
	
	mpz_t ss;
	mpz_init(ss);
	mpz_t ciphertext1;
	mpz_init(ciphertext1);
	
	printf("DGHV_ecrypt\n");
	start = clock();
	for(i=0;i<10;i++) DGHV_encrypt(ciphertext1,1,paraSet.tau,paraSet.Rho);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; duration=duration/10;
	printf("DGHV_encrypt time: %f seconds\n",duration);

	printf("DGHV_decrypt...\n");
	start = clock();
	for(i=0;i<10;i++) DGHV_decrypt(ciphertext1);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; duration=duration/10;
	printf("DGHV_deccrypt time: %f seconds\n",duration);

	
	printf("Homomorphic multiplication start.\n");
	start = clock();
	for(i=0;i<10;i++) DGHV_mul(ss,ciphertext1,ciphertext1);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; duration=duration/10;
	printf("Homomorphic multiplication time: %f seconds\n",duration);

	printf("Homomorphic addition start.\n");
	start = clock();
	for(i=0;i<10;i++) DGHV_add(ss,ciphertext1,ciphertext1);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; duration=duration/10;
	printf("Homomorphic addition time: %f seconds\n",duration);

	

	printf("private key p switch y[i] with %ld bits of precision after the binary point\n",para.gam+64);
	start = clock();
	expend_p2y(para.gam,para.theta,para.Theta);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
	printf("private key p switch y[i] success! The time:%f\n",duration);
	
	printf("The time of DGHV_genkey(paraSet) function runn over! Total:%f\n",total);
	
}



void expend_ciphertext(mpz_t ciphertext,unsigned long Theta,unsigned long gam,unsigned long prec,unsigned long tau,mp_bitcnt_t Rho,mp_bitcnt_t rho){ 

	int i,j;
	char *binary_ciphertext ;
	mp_exp_t exponent ;
	mpz_t s_ciphertext;
	mpf_t c,z;
	
	mpz_init(s_ciphertext);
	mpf_init(c);
	mpf_set_z(c,ciphertext);
	mpf_init(z);
	
	for(i=0;i<prec+1;i++){
		for(j=0;j<Theta;j++){
			mpz_init_set_ui(evla_matrix[j][i],0);
		}
	}
        srand(time(NULL));
	printf("start expend_ciphertext...\n");
	start = clock();
	binary_ciphertext=(char *)malloc((gam+128)*sizeof(char));
	for(i=0;i<Theta;i++){
		
		mpf_mul(z,y[i],c);
        	if(s[i]==1){
        		DGHV_encrypt(s_ciphertext,s[i],tau,rho); 
		
        	}else{
            		if(i%20==0){
               			DGHV_encrypt(s_ciphertext,s[i],tau,rho); 
                		printf("OK\n");
		
            		}else{

		        	mpz_set_ui(s_ciphertext,s[i]);
            		}
            		//DGHV_encrypt(s_ciphertext,s[i],tau,rho); 
		
        	}
		
		mpf_get_str(binary_ciphertext,&exponent,2,0,z);
		
		if(exponent>0){
			
			for(j=0;j<prec+1;j++){
				mpz_mul_ui(evla_matrix[i][j],s_ciphertext,binary_ciphertext[exponent-1+j]-'0');  
			}
		}else if(exponent<0){
			for(j=0;j<-exponent+1;j++){
				mpz_mul_ui(evla_matrix[i][j],s_ciphertext,0); 
			}
			for(j=-exponent+1;j<prec+1;j++){
				mpz_mul_ui(evla_matrix[i][j],s_ciphertext,binary_ciphertext[exponent-1+j]-'0');
			}
		}else if(exponent==0){
			if(mpf_cmp_d(z,0.5)>=0){
				mpz_mul_ui(evla_matrix[i][j],s_ciphertext,0); 
				for(j=1;j<prec+1;j++){
					mpz_mul_ui(evla_matrix[i][j],s_ciphertext,binary_ciphertext[j]-'0'); 
				}		
			}
		}	
	}
	free(binary_ciphertext);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; 
	printf("expend_ciphertext success! The time:%f\n",duration);
	for(i=0;i<Theta;i++){
		mpf_clear(y[i]);
	}
	
	mpz_clear(s_ciphertext);
	mpf_clear(z);
	mpf_clear(c);
}

void lsb_first_term_encrypt(mpz_t lsb_ciphertext,mpz_t old_ciphertext,unsigned long tau,mp_bitcnt_t Rho){
	unsigned long int m;
	printf("start lsb_first_term_encrypt (encrypt(c mod 2))...\n");
	mpz_fdiv_r_ui(old_ciphertext,old_ciphertext,2);
	m = mpz_get_ui(old_ciphertext);
	DGHV_encrypt(lsb_ciphertext,m,tau,Rho);
	printf("lsb_first_term_encrypt sucess!\n");
}

void lsb_second_term_encrypt(mpz_t lsb_ciphertext,unsigned long prec,unsigned long theta,unsigned long Theta){
	
	int i,j;
	int c,d,l,k,t,m;
	printf("start lsb_second_term_encrypt (encrypt(「∑si*zi」mod 2))...\n");
	
	for(i=prec;i>0;i--){
		c=(int)floor(log(theta+prec-i)/log(2));
		l=Theta+prec-i;
		d=i-c;
		if(d>=0) d=0;
		for(j=1,t=1;j<=(int)pow(2,c+d);j++){
			for(k=j;k<=l;k++){
				mpz_init_set_ui(hw_matrix[j][k],0);
				DGHV_mul(hw_matrix[j][k],evla_matrix[k-1][i],hw_matrix[j-1][k-1]);
				DGHV_add(hw_matrix[j][k],hw_matrix[j][k],hw_matrix[j][k-1]);
			}
			if(j==(int)pow(2,t)){
				mpz_set(evla_matrix[Theta+(prec-i)][prec-t],hw_matrix[j][l]);
				t++;
			}
			if(j>=2){
				for(m=j-1;m<=Theta;m++){
					mpz_clear(hw_matrix[j-1][m]);
				}
			}
		}
		for(m=j-1;m<=Theta;m++){
			mpz_clear(hw_matrix[j-1][m]);
		}
		for(m=0;m<l;m++) mpz_clear(evla_matrix[m][i]);
	}
	
	for(i=0;i<Theta+prec;i++){
		DGHV_add(lsb_ciphertext,lsb_ciphertext,evla_matrix[i][0]);
	}
	for(i=0;i<Theta+prec;i++){
        	mpz_clear(evla_matrix[i][0]);
        }
	printf(" lsb_second_term_encrypt sucess!\n");
}

void DGHV_fresh_noise(mpz_t new_ciphertext,mpz_t old_ciphertext,unsigned long Theta,unsigned long theta,unsigned long gam,unsigned long prec,unsigned long tau,mp_bitcnt_t Rho,mp_bitcnt_t rho){
	
	int i;
	mpz_t lsb1_ciphertext;
	mpz_t lsb2_ciphertext;
	mpz_init_set_ui(lsb1_ciphertext,0);
	mpz_init_set_ui(lsb2_ciphertext,0);
	
	
	expend_ciphertext(old_ciphertext,Theta,gam,prec,tau,Rho,rho);
	
	printf("Start fresh ciphertext...\n");
	
	start = clock();
	lsb_first_term_encrypt(lsb1_ciphertext,old_ciphertext,tau,Rho);
	printf("lsb1_ciphertext decrypt:%ld\n",DGHV_decrypt(lsb1_ciphertext));

	lsb_second_term_encrypt(lsb2_ciphertext,prec,theta,Theta);
	printf("lsb2_ciphertext decrypt:%ld\n",DGHV_decrypt(lsb2_ciphertext));
	DGHV_add(new_ciphertext,lsb1_ciphertext,lsb2_ciphertext);
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; 
	printf("Fresh ciphertext success! The time:%f\n",duration);
	
}

int main(){

	int i,j;
	int in;
	unsigned long m1,m2;
	
	m1=1;
	m2=0;
	
	mpz_t ciphertext1,ciphertext2,result1,result2,fresh_ciphertext;
	mpz_init(ciphertext1);
	mpz_init(ciphertext2);
	mpz_init(result1);
	mpz_init(result2);
	mpz_init(fresh_ciphertext);
	mpz_init(p);
	
	printf("choose paramters set:\n<L0>    0\n<L1>    1\n<L2>    2\n<L3>    3\n<L4>    4\n");
	printf("parameters\tlam\trho\teta\tgam\t\ttao\n");
	printf("L0\t\t52\t24\t1632\t2.0*10^6\t1000\n");
	printf("L1\t\t52\t24\t1728\t2.5*10^6\t1000\n");
	printf("L0\t\t52\t24\t1801\t3.0*10^6\t1000\n");
	printf("L0\t\t52\t24\t1874\t3.5*10^6\t1000\n");
	printf("L0\t\t52\t24\t1993\t4.2*10^6\t1000\n");
	
	scanf("%d",&in);
	if(in==L0)  DGHV_init(L0);
	if(in==L1)  DGHV_init(L1);
	if(in==L2)  DGHV_init(L2);
	if(in==L3)  DGHV_init(L3);
	if(in==L4)  DGHV_init(L4);

	DGHV_genkey(paraSet);
	
	printf("m1=%lu\n",m1);
	printf("m2=%lu\n",m2);
	
	start = clock();
	DGHV_encrypt(ciphertext1,m1,paraSet.tau,paraSet.Rho);
	//gmp_printf("ciphertext1:%Zd\n",ciphertext1);
	DGHV_encrypt(ciphertext2,m2,paraSet.tau,paraSet.Rho);
	//gmp_printf("ciphertext2:%Zd\n",ciphertext2);
	finish = clock();  
	duration = (double)(finish - start) / CLOCKS_PER_SEC;    
    printf( "DGHV_encrypt:%.15f seconds\n", duration );    
	
	start = clock();
	printf("m1 decrypt:%ld\n",DGHV_decrypt(ciphertext1));
	printf("m2 decrypt:%ld\n",DGHV_decrypt(ciphertext2));
	finish = clock();  
	duration = (double)(finish - start) / CLOCKS_PER_SEC;    
    printf( "DGHV_decrypt:%.15f seconds\n", duration ); 
	
	printf("m1+m2=1+0=1\n");
	//gmp_printf("ciphertext1:%Zd\n",ciphertext1);
	//gmp_printf("ciphertext2:%Zd\n",ciphertext2);
	DGHV_add(result1,ciphertext1,ciphertext2);
	//gmp_printf("ciphertext1+ciphertext2=%Zd\n",result1);
	printf("result1 decrypt:%ld\n\n",DGHV_decrypt(result1));
	
	printf("m1*m2=1*0=0\n");
	//gmp_printf("ciphertext1:%Zd\n",ciphertext1);
	//gmp_printf("ciphertext2:%Zd\n",ciphertext2);
	DGHV_mul(result2,ciphertext1,ciphertext2);
	//gmp_printf("ciphertext1*ciphertext2=%Zd\n",result2);
	printf("result2 decrypt:%ld\n\n",DGHV_decrypt(result2));
	
	DGHV_fresh_noise(fresh_ciphertext,result2,paraSet.Theta,paraSet.theta,paraSet.gam,paraSet.n,paraSet.tau,paraSet.Rho,paraSet.rho); 
	printf("Fresh_ciphertext decrypt:%ld\n\n",DGHV_decrypt(fresh_ciphertext));
	
	system("pause");
	return 0;
}























