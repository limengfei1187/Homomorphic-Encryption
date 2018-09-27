/*************************************************************************
    > File Name: CMNT.c
    > Author: Mengfei Li
    > Mail: 1187337437@qq.com 
    > Created Time: 2018-09-27 Thursday 11:55:50
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
    unsigned long eta;     //私钥长度
    unsigned long gam;     //公钥长度
    unsigned long Theta;   //集合维度
    unsigned long theta;   //稀疏子集长度 默认长度为15
    unsigned long tau;     //公钥个数
    unsigned long prec;       //精度 默认等于4
    unsigned long beta;    //段共要的新参数，可以降低公钥大小
    unsigned long kappa;    //kappa = gam +2
}parameters;


typedef struct publickeys{
    mpz_t x0;                        //X0
    mpz_t pk_vector1[100];         //Xi0
    mpz_t pk_vector2[100];         //Xj1
    mpz_t s0_vector[1000];         //S(0)向量
    mpz_t s1_vector[1000];         //S(1)向量
    mpz_t s_fills[1000];
    mpf_t y[10000];
    unsigned long int seed;       //随机生成Uij的种子
    int s0_group_cnt;
    int s1_group_cnt;
    int every_group_length;
    int last_group_length;
    int fill_cnt;

}publickeys;

typedef struct privatekeys{
    mpz_t p;
    int s0[30][30];
    int s1[30][30];
    int fill_s[1000];
    int s0_group_cnt;
    int s1_group_cnt;
    int every_group_length;
    int last_group_length;
    int fill_cnt;

}privatekey;

int rand_integer_matrix[100][100];
double start,finish,duration;
privatekey sk;
publickeys pk_set;
parameters para_set;
mpz_t z_base;
mpf_t f_base;
mpz_t evla_matrix[10000][20];
mpz_t hw_matrix[10][10000];
gmp_randstate_t rstate;
unsigned long seed;
struct timeval t_val;

void CMNT_init(const int para_type){
	
  if(para_type == L0){
        para_set.lam = 52;
        para_set.Rho = 24;
        para_set.rho = 24;
        para_set.eta = 1632;
        para_set.gam = 2000000;
        para_set.Theta = 500;
        para_set.theta = 15;
        para_set.prec = 4;
        para_set.beta = 32;
        para_set.kappa = 2500024;
    }
  
    if(para_type == L1){
        para_set.lam = 52;
        para_set.Rho = 24;
        para_set.rho = 24;
        para_set.eta = 1728;
        para_set.gam = 2500000;
        para_set.Theta = 500;
        para_set.theta = 15;
        para_set.prec = 4;
        para_set.beta = 32;
        para_set.kappa = 2000030;
    }
    if(para_type ==  L2){
        para_set.lam = 52;
        para_set.Rho = 24;
        para_set.rho = 24;
        para_set.eta = 1801;
        para_set.gam = 3000000;
        para_set.Theta = 500;
        para_set.theta = 15;
        para_set.prec = 4;
        para_set.beta = 32;
        para_set.kappa = 3000032;

    }
    if(para_type ==  L3){
        para_set.lam = 52;
        para_set.Rho = 24;
        para_set.rho = 24;
        para_set.eta = 1874;
        para_set.gam = 3500000;
        para_set.Theta = 500;
        para_set.theta = 15;
        para_set.prec = 4;
        para_set.beta = 32;
        para_set.kappa = 3500040;

    }
    if(para_type ==  L4){
        para_set.lam = 52;
        para_set.Rho = 24;
        para_set.rho = 24;
        para_set.eta = 1993;
        para_set.gam = 4200000;
        para_set.Theta = 500;
        para_set.theta = 15;
        para_set.prec = 4;
        para_set.beta = 32;
        para_set.kappa = 4200048;

    }
	int i,j;
	mpf_set_default_prec(para_set.kappa);
	mpf_init_set_ui(f_base,2);
	mpz_init_set_ui(z_base,2);
	for(i=0;i<=para_set.Theta;i++){
		mpz_init_set_ui(hw_matrix[0][i],1);
		mpf_init(pk_set.y[i]);
	}
	

    gettimeofday(&t_val, NULL);
    seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate,seed);

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

void gen_s(){
    int i,j;
    int t,b,sub;
    int add_cnt=0;
    int reduce_cnt=0;
    unsigned long int seed;

    printf("Theta:%ld\nbeta:%ld\n",para_set.Theta,para_set.beta);
    t = (int)floor(sqrt(para_set.theta));
    printf("theta:%ld\nsqrt(theta):%d\n",para_set.theta,t);
    if(t*1.0==sqrt(para_set.theta)){
        sk.s0_group_cnt = t;
        sk.s1_group_cnt = t;
    }
    if(para_set.theta<=t*(t+1)){
        sk.s0_group_cnt = t+1;
        sk.s1_group_cnt = t;

    }
    if(para_set.theta>t*(t+1)&&para_set.theta<(t+1)*(t+1)){
        sk.s0_group_cnt = t+1;
        sk.s1_group_cnt = t+1;
    }
    printf("s0_group_cnt:%d\ns1_group_cnt:%d\n",sk.s0_group_cnt,sk.s1_group_cnt);

    b=(int)round(sqrt(para_set.Theta/(sk.s0_group_cnt*sk.s1_group_cnt)*1.0));

    printf("b=sqrt(B):%d\n",b);

    sub=para_set.Theta-sk.s0_group_cnt*sk.s1_group_cnt*b*b; 
    if(sub==0){
        sk.every_group_length = b;
        sk.last_group_length = b;
        sk.fill_cnt = 0;
    }else if(sub < 0){
        sk.every_group_length = b;
        reduce_cnt = (int)ceil(sub*(-1)/(sk.s0_group_cnt*b*1.0));
       
        sk.last_group_length=b-reduce_cnt;
        sk.fill_cnt = reduce_cnt*sk.s0_group_cnt*b+sub;
    }else if(sub > 0){
        sk.every_group_length = b;
        add_cnt = (int)floor(sub/(sk.s0_group_cnt*b*1.0))*(-1);
        sk.last_group_length=b+add_cnt;
        sk.fill_cnt = sub+add_cnt*sk.s0_group_cnt*b;
    }
    printf("every_group_length:%d\n",sk.every_group_length);
    printf("last_group_length:%d\n",sk.last_group_length);
    printf("add_cnt:%d\n",add_cnt);
    printf("reduce_cnt:%d\n",reduce_cnt);
    printf("fill_cnt:%d\n",sk.fill_cnt);

    int length,tag;
    seed = time(NULL);
    srand(seed);
    for(i=0;i<sk.s1_group_cnt;i++){
        printf("sk.s1:[%d]  ",i);
        if(i!=sk.s1_group_cnt-1){
            length = sk.every_group_length;
        }else{
            length = sk.last_group_length;

        }
        if(i==sk.s1_group_cnt-1){
            tag=0;
            sk.s1[i][tag]=1;
        }else{

            tag = rand()%length;
            sk.s1[i][tag]=1;
        }
        for(j=0;j<length;j++){
            if(j!=tag){
                sk.s1[i][j]=0;
            }
            printf("%d  ",sk.s1[i][j]);

        }
        printf("\n");
    }


    for(i=0;i<sk.s0_group_cnt;i++){
        printf("sk.s0:[%d]  ",i);
        length = sk.every_group_length;
        tag = rand()%length;
        sk.s0[i][tag]=1;
        for(j=0;j<length;j++){
            if(j!=tag){
                sk.s0[i][j]=0;
            }
            printf("%d  ",sk.s0[i][j]);
        }
        printf("\n");
    }

    printf("sk.fill:");
    for(i=0;i<sk.fill_cnt;i++){
        sk.fill_s[i]=0;
        printf("%d ",sk.fill_s[i]);
    }
    printf("\n");
    }

void gen_privatekey(){

    mpz_t rand_num;
    mpz_init(sk.p);
    mpz_init(rand_num);
    mpz_rrandomb(rand_num,rstate,para_set.eta);
    mpz_nextprime(sk.p,rand_num);
    gen_s();
    mpz_clear(rand_num);

}

void expend_p2y(){

    int i,j,k,l,t;
	int length,cnt;
    unsigned long int seed,random,last_group_length;

    mpz_t yy[15];
    mpz_t quotient, remainder,result;
    mpz_t xp;
    mpz_t rand_num,ui;
    mpf_t numerator,denominator;

    mpz_init(quotient);
    mpz_init(remainder);
    mpz_init(result);
    mpz_init(xp);
    mpz_init(rand_num);
    mpz_init(ui);
    mpf_init(numerator);
    mpf_init(denominator);

    mpz_ui_pow_ui(xp,2,para_set.gam+2);
    div_trun(xp,xp,sk.p);
    mpz_fdiv_qr_ui(quotient,remainder,xp,para_set.theta);
    for(i=0; i<para_set.theta; i++){
        mpz_init_set_ui(yy[i],0);
        mpz_add(yy[i], yy[i], quotient);
    }
    mpz_add(yy[para_set.theta-1], yy[para_set.theta-1], remainder);

    for(i=0; i<para_set.theta; i++){
        mpz_fdiv_qr_ui(quotient,remainder,yy[i],rand()%para_set.rho+1);
        mpz_add(result,quotient,remainder);
        mpz_sub(yy[i],yy[i],result);

        mpz_fdiv_qr_ui(quotient,remainder,result,para_set.theta);
        mpz_add(result,remainder,quotient);

        for(j=0;j<para_set.theta;j++){
            if(j==i){
                mpz_add(yy[j],yy[j],result);

            }else{
                mpz_add(yy[j],yy[j],quotient);
            }
        }
    }
	
    mpz_ui_pow_ui(ui,2,para_set.gam+3);
    mpf_pow_ui(denominator,f_base,para_set.gam+2);
    
    cnt=t=0;
    for(i=0;i<sk.s0_group_cnt;i++){
        for(j=0;j<sk.s1_group_cnt;j++){
         if(j==sk.s1_group_cnt-1){
             length=sk.last_group_length;
         }else{
            length=sk.every_group_length;
        }
        for(k=0;k<sk.every_group_length;k++){
            for(l=0;l<length;l++){
                if(t<para_set.theta){
                    if(sk.s0[i][k]*sk.s1[j][l]==1){
                        mpf_set_z(numerator,yy[t]);
                        t++;
                    }else{
                        mpz_urandomm(rand_num,rstate,ui);
                        mpf_set_z(numerator,rand_num);
                    }
                }else{
                    mpz_urandomm(rand_num,rstate,ui);
                    mpf_set_z(numerator,rand_num);
                }
                mpf_div(pk_set.y[cnt],numerator,denominator);
		cnt++;
            }
        }
    }
}
for(i=0;i<sk.fill_cnt;i++){
    mpz_urandomm(rand_num,rstate,ui);
    mpf_set_z(numerator,rand_num);
    mpf_div(pk_set.y[cnt],numerator,denominator);
	cnt++;
}

for(i=0;i<para_set.theta;i++) mpz_clear(yy[i]);

    mpz_clear(rand_num);
    mpz_clear(ui);
    mpz_clear(quotient);
    mpz_clear(remainder);
    mpz_clear(result);
    mpz_clear(xp);
    mpf_clear(numerator);
    mpf_clear(denominator);

}

void gen_pk_set(){
    int i;

    mpz_t q,r,upper_b;
    mpz_init(q);
    mpz_init(r);
    mpz_init(upper_b);

    gettimeofday(&t_val, NULL);
    seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;

    mpz_pow_ui(upper_b,z_base,para_set.gam);
    div_trun(upper_b,upper_b,sk.p);
    mpz_urandomm(q,rstate,upper_b);
    mpz_mul(pk_set.x0,sk.p,q);
    mpz_set(upper_b,q);
    for(i=0;i<=para_set.beta;i++){
        mpz_init(pk_set.pk_vector1[i]);
        mpz_init(pk_set.pk_vector2[i]);
        mpz_urandomm(q,rstate,upper_b);
        mpz_rrandomb(r,rstate,para_set.Rho);
        mpz_mul(pk_set.pk_vector1[i],sk.p,q);
        mpz_add(pk_set.pk_vector1[i],pk_set.pk_vector1[i],r);
        mpz_urandomm(q,rstate,upper_b);
        mpz_rrandomb(r,rstate,para_set.Rho);
        mpz_mul(pk_set.pk_vector2[i],sk.p,q);
        mpz_add(pk_set.pk_vector2[i],pk_set.pk_vector1[i],r);

    }
    for(i=0;i<sk.s0_group_cnt*sk.every_group_length;i++){
        mpz_rrandomb(r,rstate,para_set.rho);
        mpz_urandomm(q,rstate,upper_b);
        mpz_mul_ui(r,r,2);
        mpz_add_ui(r,r,sk.s0[i/sk.every_group_length][i%sk.every_group_length]);
        mpz_mul(pk_set.s0_vector[i],sk.p,q);
        mpz_add(pk_set.s0_vector[i],pk_set.s0_vector[i],r);
    }

    for(i=0;i<(sk.s1_group_cnt-1)*sk.every_group_length;i++){
        mpz_rrandomb(r,rstate,para_set.rho);
        mpz_urandomm(q,rstate,upper_b);
        mpz_mul_ui(r,r,2);
        mpz_add_ui(r,r,sk.s1[i/sk.every_group_length][i%sk.every_group_length]);
        mpz_mul(pk_set.s1_vector[i],sk.p,q);
        mpz_add(pk_set.s1_vector[i],pk_set.s1_vector[i],r);
    }
    int t=(sk.s1_group_cnt-1)*sk.every_group_length;
    for(i=0;i<sk.last_group_length;i++){
        mpz_rrandomb(r,rstate,para_set.rho);
        mpz_urandomm(q,rstate,upper_b);
        mpz_mul_ui(r,r,2);
        mpz_add_ui(r,r,sk.s1[pk_set.s1_group_cnt-1][i]);
        mpz_mul(pk_set.s1_vector[t+i],sk.p,q);
        mpz_add(pk_set.s1_vector[t+i],pk_set.s1_vector[t+i],r);
    }
    t=(sk.s1_group_cnt-1)*sk.every_group_length+sk.last_group_length;
    for(i=0;i<sk.fill_cnt;i++){
        mpz_rrandomb(r,rstate,para_set.rho);
        mpz_urandomm(q,rstate,upper_b);
        mpz_mul_ui(r,r,2);
        mpz_add_ui(r,r,sk.fill_s[i]);
        mpz_mul(pk_set.s1_vector[t+i],sk.p,q);
        mpz_add(pk_set.s1_vector[t+i],pk_set.s1_vector[t+i],r);
    }
    pk_set.s0_group_cnt=sk.s0_group_cnt;
    pk_set.s1_group_cnt=sk.s1_group_cnt;
    pk_set.every_group_length=sk.every_group_length;
    pk_set.last_group_length=sk.last_group_length;
    pk_set.fill_cnt=sk.fill_cnt;
    pk_set.seed=seed;

    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(upper_b);
}

void save(const char* type){
    int i;
    FILE *out;
    if(strcmp(type,"publickey")==0){
        if((out = fopen("publickey","wt"))== NULL){
            fprintf(stderr,"Cannot open publickey file\n");
        }
        gmp_fprintf(out,"**********|pk_vector1|**********\n");
        for(i=0;i<para_set.beta;i++){
            gmp_fprintf(out,"%Zd\n",pk_set.pk_vector1[i]);
        }
        gmp_fprintf(out,"**********|pk_vector2|**********\n");
        for(i=0;i<para_set.beta;i++){
            gmp_fprintf(out,"%Zd\n",pk_set.pk_vector2[i]);
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

void CMNT_encrypt(mpz_t ciphertext, unsigned long plaintext){

    int i,j,r;
    mpz_t rand_num,product;

    mpz_init(rand_num);
    mpz_init(product);
    gettimeofday(&t_val, NULL);
    seed = t_val.tv_sec*1000*1000 + t_val.tv_usec;
    srand(seed);
    for(i=0;i<3;i++){

        j = rand()%para_set.beta;
        r = rand()%para_set.beta;
        mpz_mul(product,pk_set.pk_vector1[j],pk_set.pk_vector2[r]);
        mpz_add(ciphertext,ciphertext,product);
    }

    mpz_mul_ui(ciphertext,ciphertext,2);
    mpz_mod(ciphertext,ciphertext,pk_set.x0);
    mpz_rrandomb(rand_num,rstate,para_set.Rho);
    mpz_mul_ui(rand_num,rand_num,2);
    mpz_add_ui(rand_num,rand_num,plaintext);
    mpz_add(ciphertext,ciphertext,rand_num);
    mpz_clear(rand_num);
    mpz_clear(product);
}

unsigned long int CMNT_decrypt(mpz_t ciphertext){
    mpz_t plaintext;

    mpz_init(plaintext);
    mpz_mod(plaintext,ciphertext,sk.p);
    mpz_mod_ui(plaintext,plaintext,2);
    return mpz_get_ui(plaintext);
}

void CMNT_add(mpz_t sum,mpz_t c1,mpz_t c2){
    mpz_add(sum,c1,c2);
    mpz_mod(sum,sum,pk_set.x0);
}

void CMNT_mul(mpz_t product,mpz_t c1,mpz_t c2){
    mpz_mul(product,c1,c2);
    mpz_mod(product,product,pk_set.x0);
}

void CMNT_genkey(){

    int i;
    double total=0;
    printf("the time of generate key function running...\n");

    
    start = clock();
    gen_privatekey();
    finish = clock();
    printf("generate private key:sk.p:");
    gmp_printf("%Zd\n",sk.p);
    duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
    printf("generate success! The time:%f seconds\n",duration);
	
	printf("Save private key...\n");
    save("privatekey");
    printf("Save public key end\n");

    printf("generate pk set...\n");
    start = clock();
    gen_pk_set();
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
    printf("generate pk set success,total:%ld, the time: %f seconds\n",2*para_set.beta,duration);
	
    printf("Save public key...\n");
    save("publickey");
	printf("Save public key end\n");
	
    mpz_t ciphertext1;
    mpz_t ss;
    mpz_init(ss);
    mpz_init(ciphertext1);
    CMNT_encrypt(ciphertext1,1);
    printf("Homomorphic multiplication start.\n");
    start = clock();
    for(i=0;i<10;i++) CMNT_mul(ss,ciphertext1,ciphertext1);
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC; duration = duration/10;
    printf("Homomorphic multiplication time: %f seconds\n",duration);
	
    printf("Homomorphic addition start.\n");
    start = clock();
    for(i=0;i<10;i++) CMNT_add(ss,ciphertext1,ciphertext1);
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC; duration = duration/10;
    printf("Homomorphic addition time: %f seconds\n",duration);


    printf("Private key p switch y[i] with %ld bits of precision after the binary point\n",para_set.kappa);
    start = clock();
    expend_p2y();
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC; total+=duration;
    printf("Private key p switch y[i] success! The time:%f\n",duration);
    printf("The time of CMNT_genkey() function runn over! Total:%f\n",total);

}


void get_binary_str(mpf_t z,mpz_t s,int i){
	
    int j;
    char *binary_ciphertext ;	
    binary_ciphertext=(char *)malloc((para_set.kappa+para_set.eta)*sizeof(char));
    mp_exp_t exponent ;
    mpf_get_str(binary_ciphertext,&exponent,2,0,z);
    if(exponent>0){
        for(j=0;j<para_set.prec+1;j++){
            mpz_mul_ui(evla_matrix[i][j],s,binary_ciphertext[exponent-1+j]-'0');
        }
    }else if(exponent<0){
        for(j=0;j<-exponent+1;j++){
            mpz_mul_ui(evla_matrix[i][j],s,0);
        }
        for(j=-exponent+1;j<para_set.prec+1;j++){
            mpz_mul_ui(evla_matrix[i][j],s,binary_ciphertext[exponent-1+j]-'0');
        }
    }else if(exponent==0){
        if(mpf_cmp_d(z,0.5)>=0){
            mpz_mul_ui(evla_matrix[i][0],s,0);
            for(j=1;j<para_set.prec+1;j++){
                mpz_mul_ui(evla_matrix[i][j],s,binary_ciphertext[j]-'0');
            }
		}
    }
    free(binary_ciphertext);
}

void expend_ciphertext(mpz_t ciphertext){

    int i,j;
    int g,h,k,l,length,t;
    mpf_t c,z;
    mpz_t s;

    mpf_init(c);
    mpf_set_z(c,ciphertext);
    mpf_init(z);
    mpz_init(s);

	for(i=0;i<8;i++){mpz_init_set_ui(hw_matrix[i][0],0);}
	for(i=0;i<para_set.prec+1;i++){
		for(j=0;j<para_set.Theta;j++){
			mpz_init(evla_matrix[j][i]);
		}
	}
    printf("start expend_ciphertext...\n");
    start = clock();

    t=0;
    i=0;
    for(g=0;g<pk_set.s0_group_cnt;g++){
        for(h=0;h<pk_set.s1_group_cnt;h++){
            if(h==pk_set.s1_group_cnt-1){
                length=pk_set.last_group_length;
            }else{
                length=pk_set.every_group_length;
            }
            for(k=0;k<pk_set.every_group_length;k++){
                for(l=0;l<length;l++){
                    if(t<para_set.theta){
                        if(sk.s0[g][k]*sk.s1[h][l]==1){
                            mpz_mul(s,pk_set.s0_vector[g*pk_set.every_group_length+k],pk_set.s1_vector[h*pk_set.every_group_length+l]);
                            mpz_mod(s,s,pk_set.x0);
                        
                        }else{
                            if(i%20==0){
                                mpz_mul(s,pk_set.s0_vector[g*pk_set.every_group_length+k],pk_set.s1_vector[h*pk_set.every_group_length+l]);
                                mpz_mod(s,s,pk_set.x0);
                                printf("OK\n");
                        
                            }else{
                                mpz_set_ui(s,sk.s0[g][k]*sk.s1[h][l]);
                            }

                        }
                        //mpz_mul(s,pk_set.s0_vector[g*pk_set.every_group_length+k],pk_set.s1_vector[h*pk_set.every_group_length+l]);
                        //mpz_mod(s,s,pk_set.x0);
                        //mpz_mul_ui(s,s,sk.s0[g][k]);
                        //mpz_mul_ui(s,s,sk.s1[h][l]);
                        
                    }else{
                        mpz_set_ui(s,sk.fill_s[0]);
                        //mpz_set_ui(s,0);
                    } 
					mpf_mul(z,pk_set.y[i],c);
					get_binary_str(z,s,i);
                    i++;
                }
            }
            t++;
        }
    }
    

    for(j=0;j<pk_set.fill_cnt;j++){
        mpz_set_ui(s,sk.fill_s[j]);
        //mpz_set_ui(s,0);
        mpf_mul(z,pk_set.y[i],c);
        get_binary_str(z,s,i);
        i++;
    }

    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC;
    printf("expend_ciphertext success! The time:%f\n",duration);

	for(i=0;i<para_set.Theta;i++){
		mpf_clear(pk_set.y[i]);
	}
    mpf_clear(z);
    mpf_clear(c);

    }

void lsb_first_term_encrypt(mpz_t lsb_ciphertext,mpz_t old_ciphertext){
    unsigned long int m;
    printf("start lsb_first_term_encrypt (encrypt(c mod 2))...\n");
    mpz_fdiv_r_ui(old_ciphertext,old_ciphertext,2);
    m = mpz_get_ui(old_ciphertext);
    CMNT_encrypt(lsb_ciphertext,m);
    printf("lsb_first_term_encrypt sucess!\n");
}

void lsb_second_term_encrypt(mpz_t lsb_ciphertext){

    int i,j;
    int c,d,l,k,t,m;
    printf("start lsb_second_term_encrypt (encrypt(「∑si*zi」mod 2))...\n");

    for(i=para_set.prec;i>0;i--){
        c=(int)floor(log(para_set.theta+para_set.prec-i)/log(2));
        l=para_set.Theta+para_set.prec-i;
        d=i-c;
        if(d>=0) d=0;
        for(j=1,t=1;j<=(int)pow(2,c+d);j++){
            for(k=j;k<=l;k++){
                mpz_init_set_ui(hw_matrix[j][k],0);
                CMNT_mul(hw_matrix[j][k],evla_matrix[k-1][i],hw_matrix[j-1][k-1]);
                CMNT_add(hw_matrix[j][k],hw_matrix[j][k],hw_matrix[j][k-1]);
            }
            if(j==(int)pow(2,t)){
                mpz_set(evla_matrix[para_set.Theta+(para_set.prec-i)][para_set.prec-t],hw_matrix[j][l]);
                t++;
            }
            if(j>=2){
                for(m=j-1;m<=para_set.Theta;m++){
                    mpz_clear(hw_matrix[j-1][m]);
                }
            }
        }
        for(m=j-1;m<=para_set.Theta;m++){
            mpz_clear(hw_matrix[j-1][m]);
        }
	for(m=0;m<l;m++) mpz_clear(evla_matrix[m][i]);
    }

    for(i=0;i<para_set.Theta+para_set.prec;i++){
        CMNT_add(lsb_ciphertext,lsb_ciphertext,evla_matrix[i][0]);
    }
    for(i=0;i<para_set.Theta+para_set.prec;i++){
        mpz_clear(evla_matrix[i][0]);
    }
    printf("lsb_second_term_encrypt sucess!\n");
}

void CMNT_fresh_noise(mpz_t new_ciphertext,mpz_t old_ciphertext){

    mpz_t lsb1_ciphertext;
    mpz_t lsb2_ciphertext;
    mpz_init_set_ui(lsb1_ciphertext,0);
    mpz_init_set_ui(lsb2_ciphertext,0);
    expend_ciphertext(old_ciphertext);

    printf("start fresh ciphertext...\n");
    printf("evalute lsb of c in ciphertext...\n");
    start = clock();
    lsb_first_term_encrypt(lsb1_ciphertext,old_ciphertext);
    printf("lsb of c decrypt in ciphertext:%ld\n",CMNT_decrypt(lsb1_ciphertext));
    printf("evalute lsb of c/p in ciphertext...\n");
    lsb_second_term_encrypt(lsb2_ciphertext);
    printf("lsb of c/p decrypt:%ld\n",CMNT_decrypt(lsb2_ciphertext));
    printf("evalute lsb(c)+lsb(c/p) in ciphertext...\n");
    CMNT_add(new_ciphertext,lsb1_ciphertext,lsb2_ciphertext);
    finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC;
    printf("fresh ciphertext success! The time:%f\n",duration);

}

int main(){

    int i,j;
    int in;
    unsigned long m1,m2;
    mpz_t ciphertext1,ciphertext2,result1,result2,fresh_ciphertext;
    m1=1;
    m2=0;
    mpz_init_set_ui(ciphertext1,0);
    mpz_init_set_ui(ciphertext2,0);
    mpz_init(result1);
    mpz_init(result2);
    mpz_init(fresh_ciphertext);
	
	printf("choose paramters set:\n<L0>    0\n<L1>    1\n<L2>    2\n<L3>    3\n<L4>    4\n");
	printf("parameters\tlam\trho\teta\tgam\t\ttao\n");
	printf("L0\t\t52\t24\t1632\t2.0*10^6\t1000\n");
	printf("L1\t\t52\t24\t1728\t2.5*10^6\t1000\n");
	printf("L0\t\t52\t24\t1801\t3.0*10^6\t1000\n");
	printf("L0\t\t52\t24\t1874\t3.5*10^6\t1000\n");
	printf("L0\t\t52\t24\t1993\t4.2*10^6\t1000\n");
	
	scanf("%d",&in);
	if(in==L0)  CMNT_init(L0);
	if(in==L1)  CMNT_init(L1);
	if(in==L2)  CMNT_init(L2);
	if(in==L3)  CMNT_init(L3);
	if(in==L4)  CMNT_init(L4);
	
    CMNT_genkey();

    printf("m1=%lu\n",m1);
    printf("m2=%lu\n",m2);
    printf("CMNE_encrypt start...\n");
    start = clock();

    CMNT_encrypt(ciphertext1,m1);
    //gmp_printf("ciphertext1:%Zd\n",ciphertext1);
    CMNT_encrypt(ciphertext2,m2);
    //gmp_printf("ciphertext2:%Zd\n",ciphertext2);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "CMNT_encrypt over:%f seconds\n", duration );

    printf("CMNT_decrypt start...\n");
    start = clock();

    printf("m1 decrypt:%lu\n",CMNT_decrypt(ciphertext1));
    printf("m2 decrypt:%lu\n",CMNT_decrypt(ciphertext2));

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "CMNT_decrypt over:%f seconds\n", duration );

    printf("m1+m2=1+0=1\n");
    printf("add evalute start...\n");
    //gmp_printf("ciphertext1:%Zd\n",ciphertext1);
    //gmp_printf("ciphertext2:%Zd\n",ciphertext2);
    CMNT_add(result1,ciphertext1,ciphertext2);
    //gmp_printf("ciphertext1+ciphertext2=%Zd\n",result1);
    printf("add evalute over,result1=CMNT_encrypt(m1)+CMNT_encrypt(m2)\n");
    printf("result1 decrypt:%lu\n\n",CMNT_decrypt(result1));

    printf("m1*m2=1*0=0\n");
    printf("multi evalute start...\n");
    //gmp_printf("ciphertext1:%Zd\n",ciphertext1);
    //gmp_printf("ciphertext2:%Zd\n",ciphertext2);
    CMNT_mul(result2,ciphertext1,ciphertext2);
    //gmp_printf("ciphertext1*ciphertext2=%Zd\n",result2);
    printf("multi evalute over,result2=CMNT_encrypt(m1)+CMNT_encrypt(m2)\n");
    printf("result2 decrypt:%lu\n\n",CMNT_decrypt(result2));

    printf("fresh noise of result2(c=result2)\n");
    CMNT_fresh_noise(fresh_ciphertext,result2);
    printf("Fresh_ciphertext decrypt:%lu\n\n",CMNT_decrypt(fresh_ciphertext));
    system("pause");
    return 0;

}

