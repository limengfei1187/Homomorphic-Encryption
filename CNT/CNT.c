/*************************************************************************
    > File Name: CMNT.c
    > Author: Mengfei Li
    > Mail: 1187337437@qq.com 
    > Created Time: 2018-09-27 Thursday 11:54:31
 ************************************************************************/

#include <stdlib.h>
#include <stdio.h>
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
    unsigned long eta;     //私钥长度
    unsigned long gam;     //公钥长度
    unsigned long Theta;   //集合维度
    unsigned long theta;   //稀疏子集长度 默认长度为15
    unsigned long tau;     //公钥个数
    unsigned long n;       //精度 默认等于4
    unsigned long beta;    //段共要的新参数，可以降低公钥大小
    unsigned long kappa;    //kappa = gam +2
}parameters;

typedef struct publickeys{
    mpz_t x0;  
    mpz_t rx0;
    mpz_t correction[8000]; 
    mpf_t y[8000];
    mpz_t sigma[3000][8000];
    mpf_t s_expand[3000][8000];
    unsigned long int seed;       
}publickeys;

typedef struct privatekeys{
    mpz_t p;
    mpz_t reduce_p;
    int s[8000];

}privatekey;


double start,finish,duration;
unsigned long seed;
struct timeval t_val;
privatekey sk;
publickeys pk;
parameters para;
mpz_t z_base;
mpf_t f_base;
gmp_randstate_t rstate;
int c_expand[8000][3000];

void CNT_init(const int para_type){

    if(para_type == L0){
        para.lam = 52;
        para.Rho = 24;
        para.rho = 24;
        para.eta = 1632;
        para.gam = 2000000;	
        para.Theta = 500;
        para.theta = 15;
        para.n = 4;
        para.tau = 1000;
        para.kappa = 2000052;
    }
	if(para_type == L1){
        para.lam = 52;
        para.Rho = 24;
        para.rho = 24;
        para.eta = 1728;
        para.gam = 2500000;	
        para.Theta = 500;
        para.theta = 15;
        para.n = 4;
        para.tau = 1000;
        para.kappa = 2001728;
    }
    if(para_type ==  L2){
        para.lam = 52;
        para.Rho = 24;
        para.rho = 24;
        para.eta = 1801;
        para.gam = 3000000;	
        para.Theta = 500;
        para.theta = 15;
        para.n = 4;
        para.tau = 1000;
        para.kappa = 3001865;

    }
    if(para_type ==  L3){
        para.lam = 52;
        para.Rho = 24;
        para.rho = 24;
        para.eta = 1874;
        para.gam = 3500000;	
        para.Theta = 500;
        para.theta = 15;
        para.n = 4;
        para.tau = 1000;
        para.kappa = 4002073;

    }
    if(para_type ==  L4){
        para.lam = 52;
        para.Rho = 24;
        para.rho = 24;
        para.eta = 1993;
        para.gam = 4200000;	
        para.Theta = 500;
        para.theta = 15;
        para.n = 4;
        para.tau = 1000;
        para.kappa = 5502489;

    }
    mpf_set_default_prec(para.kappa);
    mpf_init_set_ui(f_base,2);
    mpz_init_set_ui(z_base,2);

    gettimeofday(&t_val, NULL);
    seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate,seed);

}


void mpf_round2_mpz(mpz_t rop,mpf_t op){
    int i,length,diff;
    mpf_t a,b;
    char *str;
    mpf_init(a);
    mpf_init(b);
    str=(char *)malloc((para.kappa)*sizeof(char));
    mp_exp_t exponent;

    mpf_trunc(a,op);
    mpf_sub(b,op,a);
    if(mpf_cmp_d(b,0.5)>=0){
        mpf_add_ui(a,a,1);
        mpf_get_str(str,&exponent,2,0,a);
    }else{
        mpf_get_str(str,&exponent,2,0,a);
    }
    length = strlen(str);
    
    
    for(i=0;i<exponent-length;i++){
        strcat(str,"0");
    }
    mpz_set_str(rop,str,2);
    free(str);
    mpf_clear(a);
    mpf_clear(b);
}

void mpf_div_round(mpz_t q, mpz_t  n, mpz_t d){

    mpf_t a,b,c;

    mpf_init(a);
    mpf_init(b);
    mpf_init(c);

    mpf_set_z(a,n);
    mpf_set_z(b,d);
    mpf_div(c,a,b);
    mpf_round2_mpz(q,c);

    mpf_clear(a);
    mpf_clear(b);
    mpf_clear(c);
}



void gen_privatekey(){

    mpz_t rand_num;
    mpz_init(sk.p);
    mpz_init(sk.reduce_p);
    mpz_init(rand_num);
    mpz_rrandomb(rand_num,rstate,para.eta);
    mpz_nextprime(sk.p,rand_num);

    mpz_rrandomb(rand_num,rstate,para.eta-para.Rho);
    mpz_nextprime(sk.reduce_p,rand_num);
    mpz_clear(rand_num);

}

void gen_pk_set(){

    int i;
    mpz_t delta,random;
    mpz_t bound_kai,bound_ksi,ksi;
    mpz_t kai[10000];
    mpz_t q;


    mpz_init(q);
    mpz_init(bound_kai);
    mpz_init(delta);
    mpz_init(random);
    mpz_init(bound_ksi);
    mpz_init(ksi);
    mpz_init(pk.x0);

    gettimeofday(&t_val, NULL);
    seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;
    gmp_randseed_ui(rstate,seed);

    pk.seed=seed;

    mpz_pow_ui(bound_ksi,z_base,para.lam+para.eta);
    mpz_pow_ui(bound_kai,z_base,para.gam);
    mpz_fdiv_q(bound_ksi,bound_ksi,sk.p);

    for(i=0;i<para.tau;i++){
        mpz_init_set_ui(pk.correction[i],0);
        mpz_init(kai[i]);
        mpz_urandomm(kai[i],rstate,bound_kai);
    }
    for(i=0;i<para.tau;i++){

        mpz_urandomm(ksi,rstate,bound_ksi);
        mpz_mul(ksi,ksi,sk.p);
        mpz_mod(kai[i],kai[i],sk.p);
        mpz_rrandomb(random,rstate,para.rho);

        mpz_add(pk.correction[i],pk.correction[i],kai[i]);
        mpz_add(pk.correction[i],pk.correction[i],ksi);
        mpz_sub(pk.correction[i],pk.correction[i],random);
    }

    mpz_pow_ui(q,z_base,para.gam);
    mpz_fdiv_q(q,q,sk.p);
    mpz_urandomm(q,rstate,q);
    if(mpz_odd_p(q)==0){
        mpz_add_ui(q,q,1);
    }
    mpz_mul(pk.x0,q,sk.p);
    
    mpz_clear(q);
    mpz_clear(delta);
    mpz_clear(random);
    mpz_clear(bound_ksi);
    mpz_clear(bound_kai);
    mpz_clear(ksi);
    for(i=0;i<para.tau;i++) mpz_clear(kai[i]);
} 


void expend_p2y(){

    int i,j,r;
    unsigned long int se;
    mpz_t yy[15];
    mpz_t quotient, remainder,result_t;
    mpf_t xp;
    mpz_t rand_num,ui;
    mpf_t numerator,denominator,result_f;


    mpz_init(quotient);
    mpz_init(remainder);
    mpz_init(result_t);
    mpf_init(xp);
    mpz_init(rand_num);
    mpz_init(ui);

    mpf_init(numerator);
    mpf_init(denominator);
    mpf_init(result_f);

 
    for(i=0;i<para.Theta;i++){ 
        sk.s[i]=0;
    }

    se = time(NULL);
    srand(se);
    
    for(i=0;i<para.theta;i++){
        r=rand()%para.Theta;
        if(sk.s[r]==0){
            sk.s[r]=1;
            mpz_init_set_ui(yy[i],0); 
        }else{

            i--;
        }
    }
    mpf_pow_ui(numerator,f_base,para.eta-para.Rho);    
    mpf_set_z(denominator,sk.p);
    mpf_div(result_f,numerator,denominator);
    mpf_pow_ui(xp,f_base,para.gam+4);
    mpf_mul(result_f,result_f,xp);
    mpf_round2_mpz(result_t,result_f);
    mpz_fdiv_qr_ui(quotient,remainder,result_t,para.theta);
    
    for(i=0; i<para.theta; i++){
        mpz_add(yy[i], yy[i], quotient);
    }
    mpz_add(yy[para.theta-1], yy[para.theta-1], remainder);
    for(i=0; i<para.theta; i++){

        mpz_fdiv_qr_ui(quotient,remainder,yy[i],rand()%para.lam+1);
        mpz_add(result_t,quotient,remainder);
        mpz_sub(yy[i],yy[i],result_t);

        mpz_fdiv_qr_ui(quotient,remainder,result_t,para.theta);
        mpz_add(result_t,remainder,quotient);

        for(j=0;j<para.theta;j++){
            if(j==i){
                mpz_add(yy[j],yy[j],result_t);

            }else{
                mpz_add(yy[j],yy[j],quotient);
            }
        }
    
    }
    mpz_set_ui(result_t,0);
    mpz_ui_pow_ui(ui,2,para.gam-para.eta);
    mpf_pow_ui(denominator,f_base,para.gam+4);

    for(i=0,j=0;i<para.Theta;i++){
        mpf_init(pk.y[i]);
        if(sk.s[i]==0){
            mpz_urandomm(rand_num,rstate,ui);
            mpf_set_z(numerator,rand_num);
            
        }else if(sk.s[i]==1){
            mpf_set_z(numerator,yy[j]);
            mpz_clear(yy[j]);
            j++;
        }
        mpf_div(pk.y[i],numerator,denominator);
    }    
    
    mpz_clear(rand_num);
    mpz_clear(ui);
    mpz_clear(quotient);
    mpz_clear(remainder);
    mpz_clear(result_t);
    mpf_clear(result_f);
    mpf_clear(xp);
    mpf_clear(numerator);
    mpf_clear(denominator);
    }


void BitDecomp(mpz_t z[],unsigned long length, unsigned long k){
    int i,j;
    int binary_length;
    char *binary_str;   
    binary_str = (char *)malloc((k+2)*sizeof(char));
    for(i=0;i<length;i++){
        mpz_get_str(binary_str,2,z[i]);
        mpz_clear(z[i]);
        binary_length = strlen(binary_str);
        if(binary_length<k&&binary_length>=0){
            for(j=0;j<k-binary_length;j++){
                c_expand[i][j]=0;
            }
            for(j=k-binary_length;j<k;j++){
                c_expand[i][j]=binary_str[j-k+binary_length]-'0';
            }
        }else if(binary_length>=k){
            for(j=0;j<k;j++){
                c_expand[i][j]=binary_str[j]-'0'; 
            }
        }
       
    } 
    free(binary_str);
}


void Powersof2(int s[],unsigned long length,unsigned long k){
    int i,j;
    mpf_t weight;
    mpf_init(weight);
    for(i=k-1;i>=0;i--){
        mpf_pow_ui(weight,f_base,i);
        for(j=0;j<length;j++){
            if(s[j]==1){
                mpf_init_set(pk.s_expand[k-i-1][j],weight);
            }else{
                mpf_init_set_ui(pk.s_expand[k-i-1][j],0);
            }
        }
    }
    mpf_clear(weight);
}

void switch_key_gen(){

    int i,j;

    mpz_t upper_b,rand_num;
    mpz_t q,d;
    mpf_t a,b,c;

    mpz_init(pk.rx0);
    mpz_init(upper_b);
    mpz_init(rand_num);
    mpz_init(q);
    mpz_init(d);
    mpf_init(a);
    mpf_init(b);
    mpf_init(c);

    mpf_set_z(a,sk.reduce_p);
    mpf_pow_ui(b,f_base,para.eta-para.Rho+1);
    mpf_div(c,a,b);

    Powersof2(sk.s,para.Theta,(para.eta-para.Rho+1));
    
    mpz_pow_ui(upper_b,z_base,para.gam);
    mpz_fdiv_q(upper_b,upper_b,sk.reduce_p);
    mpz_urandomm(q,rstate,upper_b);
    mpz_mul(pk.rx0,q,sk.reduce_p);
    mpz_set(upper_b,q);

    srand(time(NULL));
    for(i=0;i<(para.eta-para.Rho+1);i++){
        for(j=0;j<para.Theta;j++){
            mpz_init(pk.sigma[i][j]);
            if(j==0&&sk.s[j]==1){
                mpz_urandomm(q,rstate,upper_b);
                mpz_rrandomb(rand_num,rstate,para.rho);
                mpz_mul(pk.sigma[i][j],sk.reduce_p,q);
                mpz_add(pk.sigma[i][j],pk.sigma[i][j],rand_num);
   
            }else{
                /*if(rand()%100101==0){
                    mpz_urandomm(q,rstate,upper_b);
                    mpz_rrandomb(rand_num,rstate,para.rho);
                    mpz_mul(pk.sigma[i][j],sk.reduce_p,q);
                    mpz_add(pk.sigma[i][j],pk.sigma[i][j],rand_num);
   
                }else{*/
                    mpz_set_ui(pk.sigma[i][j],0);
           // }
            }
            
            mpf_mul(c,c,pk.s_expand[i][j]);
            mpf_round2_mpz(d,c);
            mpz_add(pk.sigma[i][j],pk.sigma[i][j],d);
        }
    }
    mpz_clear(upper_b);
    mpz_clear(rand_num);
    mpz_clear(q);
    mpz_clear(d);
    mpf_clear(a);
    mpf_clear(b);
    mpf_clear(c);
}

void key_switch(mpz_t new_ciphertext,mpz_t old_ciphertext){
    int i,j;
    mpz_t z[10000];
    mpf_t zf,zz;
    mpz_t lsb,product;
    mpz_t upper_b;
    mpz_init(lsb);
    mpz_init(product);
    mpz_init(upper_b);
    mpf_init(zf);
    mpf_init(zz);
    mpz_pow_ui(upper_b,z_base,para.eta-para.Rho+1);
    mpf_set_z(zf,old_ciphertext);
    mpz_set_ui(new_ciphertext,0);
	
    printf("Expand ciphertext start...\n");
    for(i=0;i<para.Theta;i++){
        mpz_init(z[i]);
        mpf_mul(zz,zf,pk.y[i]);
        mpf_round2_mpz(z[i],zz);
        mpz_mod(z[i],z[i],upper_b);
        //mpf_set_z(zf,old_ciphertext);
    }
    printf("BitDecomp start\n");
    BitDecomp(z,para.Theta,(para.eta-para.Rho+1));
    printf("BitDecomp end\n");
	
    for(i=0;i<(para.eta-para.Rho+1);i++){
        for(j=0;j<para.Theta;j++){
            mpz_mul_ui(product,pk.sigma[i][j],c_expand[j][i]);
            mpz_add(new_ciphertext,new_ciphertext,product);
            mpz_mul_ui(new_ciphertext,new_ciphertext,2);
            mpz_mod(new_ciphertext,new_ciphertext,pk.rx0);
            mpz_clear(pk.sigma[i][j]);
            //mpz_clear(c_expand[j][i]);
        }
    }
    //mpz_mul_ui(product,product,2);
    mpz_mod_ui(lsb,old_ciphertext,2);
    mpz_add(new_ciphertext,product,lsb);
    //mpz_mod(new_ciphertext,new_ciphertext,pk.rx0);

    mpf_clear(zf);
    mpz_clear(lsb);
    mpz_clear(product);
    mpz_clear(upper_b);
    //for(i=0;i<para.theta;i++) mpz_clear(z[i]);

}

void swap(int *a, int *b){  
    int temp;  
  
    temp = *a;  
    *a = *b;  
    *b = temp;  
  
    return ;  
}  

void quicksort(int array[], int maxlen, int begin, int end){
    int i, j;
    if(begin < end){
        i = begin + 1; 
        j = end;        
          
        while(i < j){
            if(array[i] > array[begin])  {
                swap(&array[i], &array[j]);  
                j--;
            }
            else{
                i++;  
            }
        }

        if(array[i] >= array[begin])  {
            i--;
        }

        swap(&array[begin], &array[i]);  
        quicksort(array, maxlen, begin, i);
        quicksort(array, maxlen, j, end);
    }
}


void CNT_encrypt(mpz_t ciphertext, unsigned long plaintext){

    int i,j;
    int r;
    unsigned long se;
    mpz_t rand_num,bound,xi,X;

    mpz_init(rand_num);
    mpz_init(bound);
    mpz_init(xi);
    mpz_init(X);
    mpz_set_ui(ciphertext,0);
   
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate,pk.seed);
    gettimeofday(&t_val, NULL);
    se = t_val.tv_sec*1000*1000 + t_val.tv_usec;
    srand(se);
   
    mpz_pow_ui(bound,z_base,para.gam);
    for(i=0;i<para.tau;i++){
        r=rand()%para.lam;
        mpz_urandomm(xi,rstate,bound);
        mpz_sub(X,xi,pk.correction[i]);
        mpz_mul_ui(X,X,r);
        mpz_add(ciphertext,ciphertext,X);

    }
      
    mpz_mul_ui(ciphertext,ciphertext,2);
    mpz_mod(ciphertext,ciphertext,pk.x0);
    mpz_rrandomb(rand_num,rstate,para.Rho);
    mpz_mul_ui(rand_num,rand_num,2);
    mpz_add_ui(rand_num,rand_num,plaintext);
    mpz_add(ciphertext,ciphertext,rand_num);

    mpz_clear(rand_num);
    mpz_clear(bound);
    mpz_clear(xi);
    mpz_clear(X);
}

unsigned long int CNT_decrypt(mpz_t ciphertext,int l){
    mpz_t plaintext; 
    mpz_init(plaintext);
    if(l==0){   
        mpz_mod(plaintext,ciphertext,sk.p);
    }
    if(l==1){
        mpz_mod(plaintext,ciphertext,sk.reduce_p);
    }  
    mpz_mod_ui(plaintext,plaintext,2);
    return mpz_get_ui(plaintext);

}

void CNT_add(mpz_t sum,mpz_t c1,mpz_t c2){
    mpz_add(sum,c1,c2);
    mpz_mod(sum,sum,pk.x0);
}

void CNT_mul(mpz_t product,mpz_t c1,mpz_t c2){
    mpz_mul(product,c1,c2);
    mpz_mod(product,product,pk.x0);
}


void save(const char* type){
    int i;
    FILE *out;
    if(strcmp(type,"publickey")==0){
        if((out = fopen("publickey","wt"))== NULL){
            fprintf(stderr,"Cannot open publickey file\n");
        }
        gmp_fprintf(out,"**********|pk|**********\n");
        for(i=0;i<para.tau;i++){
            gmp_fprintf(out,"%Zd\n",pk.correction[i]);
        }
        
    }
    if(strcmp(type,"privatekey")==0){
        if((out = fopen("privatekey","wt"))== NULL){
            fprintf(stderr,"Cannot open publickey file\n");
        }
        gmp_fprintf(out,"**********|privatekey.p|**********\n");
        gmp_fprintf(out,"%Zd\n",sk.p);
    }
   
    fclose(out);

}
int main(){

    int in,i;
    int m1,m2;
	double total=0;
	
    mpz_t c1,c2;
    mpz_t result1, result2;
    mpz_t new_ciphertext;

    mpz_init(c1);
    mpz_init(c2);
    mpz_init(result1);
    mpz_init(result2);
    mpz_init(new_ciphertext);

    printf("choose paramters\n");
	printf("parameters\tlam\trho\teta\tgam\t\ttao\n");
	printf("L0\t\t52\t24\t1632\t2.0*10^6\t1000\n");
	CNT_init(L0);
	
	printf("Generate private key p start\n");
    start = clock();
    gen_privatekey();
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
	printf("Generate private key end. time:%f\n",duration);
	
	printf("Save private key...");
	save("privatekey");
	printf("Save private key sucessful.\n");
	
	
	printf("Generate public key set start\n");
	start = clock();
    gen_pk_set();
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
	printf("Generate public key end. time:%f\n",duration);
	
	printf("Save public key set...");
	save("publickey");
	printf("Save public key set sucessful.\n");
    
	printf("Private key p switch y[i] with %ld bits of precision after the binary point\n",para.kappa);
    start = clock();
    expend_p2y();
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
	printf("Expend_p2y end, time:%f\n",duration);
    
    printf("Genkey over, total time:%f seconds\n",total);
	
    
    m1=0;
    m2=1;
    printf("m1=%d\n",m1);
    printf("m2=%d\n",m2);

    printf("Start encrypt m1 and m2...\n");
    start = clock();

    CNT_encrypt(c1,m1);
    CNT_encrypt(c2,m2);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("CNT_encrypt over:%f seconds\n",duration);

    printf("Start decrypt m1 and m2...\n");
    start = clock();

    printf("m1 decrypt is:%lu\n",CNT_decrypt(c1,0));
    printf("m2 decrypt is:%lu\n",CNT_decrypt(c2,0));

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("CNT_decrypt is over:%f\n",duration);
   
    printf("Start homomorphic addition.\n");
    printf("m1+m2=%d\n",m1+m2); 
    start = clock();
    for(i=0;i<10;i++) CNT_add(result1,c1,c2);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;

    printf("Result1 decrypt is%lu\n",CNT_decrypt(result1,0));
    printf("Homomorphic addition is over:%f\n",duration/10);

    printf("Start homomorphic multiplication.\n");
    printf("m1*m2=%d\n",m1*m2);
    start = clock();
    for(i=0;i<10;i++) CNT_mul(result2,c1,c2);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("result2 decrypt is%lu\n",CNT_decrypt(result2,0));
    printf("Homomorphic multiplication is over:%f seconds\n",duration/10);

    printf("Generator switch key...\n");
    start = clock();
    switch_key_gen();
    finish = clock();
    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("Generator switch key over:%f\n",duration);

    printf("Start key switch...\n");
    start = clock();

    key_switch(new_ciphertext,result2);

    finish = clock();

    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Key switch over:%f\n",duration);

    printf("New ciphertext decrypt is %lu\n",CNT_decrypt(new_ciphertext,1));
    system("pause");
}

