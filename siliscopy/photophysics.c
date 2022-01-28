#include<stdio.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<ctype.h>

char ***X;
double ***K, dt, ***Prob;
int nlams=1, lams[10], *n_states, *ntot_states, *n_fluor_states;

int ***Nfs, ***Trans, **Num_trans,silent, ***max_Nfs;
//int len_var[4][]; //X, K, Kf, Nf,
char lam_atoms[10][20][10];
int nlam_atoms[10];

/* funciton: strcmp2
 * -----------------
 * Compares part of two strings.
 *
 *     *a: String a
 *     *b: String b
 *     n1: Character number to begin comparing
 *     n2: Character number to end comparing 
 */
int strcmp2(char *a, char *b, int n1,int n2){
//Compares string a[n1:n2] with string b[0:n2-n1] 
    int i, res;
    //0 is success
    for(i=n1;i<=n2;i++){
        res=abs(a[i] - b[i-n1]);
        if (res!=0){// character does not match
            return -1;
        }
    }
    return 0;
}



int get_index(char word[6], int index){
    int cnt=0;
    while(1){
        if(strcmp(word,X[index][cnt])==0){
            return cnt;
        }
        cnt++;
    }
    return -1; //not found
}


void ini_trans(char filename[50]){
    FILE *f;
    f=fopen(filename,"r");
    if (f==NULL){
        printf("File %s does not exist\n",filename);
        exit(-1);
    }
    char line[2000],*word;
    int i,j,k,cnt,cnt2,cnt3,s1,s2,fluor,foo;
    double val;

    cnt=0;
    cnt2=-1;
    while (fgets(line,2000,f) != NULL){
        word=strtok(line," \n");
        if (word==NULL){ //empty line
            continue;
        }
        if (cnt==0){
            foo=atoi(word);
            if (nlams!=foo){
	        printf("nlams %d foo %d",nlams,foo);
                printf("Mismatch of nlams in photophyics and parameters file!\n Using nlams = %d\n",nlams);
            }
            n_states=(int *)malloc(nlams*sizeof(int));
            ntot_states=(int *)malloc(nlams*sizeof(int));
            n_fluor_states=(int *)malloc(nlams*sizeof(int));
            Nfs=(int ***)malloc(nlams*sizeof(int **));
            max_Nfs=(int ***)malloc(nlams*sizeof(int **));
            X=(char ***)malloc(nlams*sizeof(char **));
            K=(double ***)malloc(nlams*sizeof(double **));
            Trans=(int ***)malloc(nlams*sizeof(int **));
            Prob=(double ***)malloc(nlams*sizeof(double **));
            Num_trans=(int **)malloc(nlams*sizeof(int *));


            word=strtok(NULL," \n");
            dt=atof(word);
            cnt=1;
        }
        else if (strchr(line,'[')!=NULL){
            word=strtok(NULL," ]");
            foo=atoi(word);
            if (lams[cnt-1]!=foo){
                printf("Error! the order of lams in parameters file and photophysics are different\n");
                printf("%d != %d\n",lams[cnt-1],foo);
                exit(-1);
            }
            lams[cnt-1]=foo;
            cnt2=0;
            cnt++;
        }            
        else if (cnt2==0){
            n_states[cnt-2]=atoi(word);
            word=strtok(NULL," \n");
            n_fluor_states[cnt-2]=atoi(word);

            ntot_states[cnt-2]=n_states[cnt-2]+n_fluor_states[cnt-2];
            Nfs[cnt-2]=(int **)malloc(n_fluor_states[cnt-2]*sizeof(int *));
            max_Nfs[cnt-2]=(int **)malloc(n_fluor_states[cnt-2]*sizeof(int *));
            for (i=0;i<n_fluor_states[cnt-2];i++){
                Nfs[cnt-2][i]=(int *)malloc(2*sizeof(int));
                max_Nfs[cnt-2][i]=(int *)malloc(2*sizeof(int));
            }
            X[cnt-2]=(char **)malloc(ntot_states[cnt-2]*sizeof(char*));
            for (i=0;i<ntot_states[cnt-2];i++){
                X[cnt-2][i]=(char *)malloc(6);
            }
            
            K[cnt-2]=(double **)malloc(ntot_states[cnt-2]*sizeof(double *));
            for(i=0;i<ntot_states[cnt-2];i++){
                K[cnt-2][i]=(double *)malloc(ntot_states[cnt-2]*sizeof(double));
                for (j=0;j<ntot_states[cnt-2];j++){
                    K[cnt-2][i][j]=0.0f;
                }
            }
            Trans[cnt-2]=(int **)malloc(ntot_states[cnt-2]*sizeof(int *));
            for(i=0;i<ntot_states[cnt-2];i++){
                Trans[cnt-2][i]=(int *)malloc(ntot_states[cnt-2]*sizeof(int));
            }

            Prob[cnt-2]=(double **)malloc(ntot_states[cnt-2]*sizeof(double *));
            for(i=0;i<ntot_states[cnt-2];i++){
                Prob[cnt-2][i]=(double *)malloc(ntot_states[cnt-2]*sizeof(double));
            }

            Num_trans[cnt-2]=(int *)malloc(ntot_states[cnt-2]*sizeof(int));
            for(i=0;i<ntot_states[cnt-2];i++){
                Num_trans[cnt-2][i]=0;
            }
            cnt2++;
        }
        else if(cnt2==1){
            strcpy(X[cnt-2][0],word);
            cnt3=1;
            while (1){
                word=strtok(NULL," \n");
                if (word==NULL)
                    break;

                strcpy(X[cnt-2][cnt3],word);
                cnt3++;
            }
            cnt2++;
        }
        else if (cnt2>1){
            s1=get_index(word,cnt-2);
            word=strtok(NULL," \n");
            s2=get_index(word,cnt-2);
            word=strtok(NULL," \n");
            val=atof(word);
            if(val<-0.5){
	        K[cnt-2][s1][s2]=-1.0f;
            }
	    else{
                K[cnt-2][s1][s2]=val*dt;
	    }
            Trans[cnt-2][s1][Num_trans[cnt-2][s1]]=s2;
            Num_trans[cnt-2][s1]++;
            if(s1>=n_states[cnt-2]){
                Nfs[cnt-2][s1-n_states[cnt-2]][0]=0;
                max_Nfs[cnt-2][s1-n_states[cnt-2]][0]=0;
                word=strtok(NULL," \n");
                Nfs[cnt-2][s1-n_states[cnt-2]][1]=atoi(word);//lam
                max_Nfs[cnt-2][s1-n_states[cnt-2]][1]=atoi(word);//lam
            }
        }
    }
    fclose(f);
    
    printf("\n\nnlams = %d\n",nlams);
    for (i=0;i<nlams;i++){
        printf("Wavelength: %d\n",lams[i]);

        for(j=0;j<ntot_states[i];j++){
            if (Num_trans[i][j]<1)
                continue;
            Prob[i][j][0]=K[i][j][Trans[i][j][0]];
            for(k=1;k<Num_trans[i][j];k++){
                Prob[i][j][k]=Prob[i][j][k-1]+K[i][j][Trans[i][j][k]];
            }
	    printf("Prob = ");
            for(k=0;k<Num_trans[i][j];k++){
                Prob[i][j][k]/=Prob[i][j][Num_trans[i][j]-1];
                printf("%f ",Prob[i][j][k]);
            }
	    printf("\n");
        }
        printf("n_states = %d\n",n_states[i]);
        printf("n_fluor_states = %d\n",n_fluor_states[i]);
        printf("X = ");
        for (j=0;j<ntot_states[i];j++){
            printf("%s ",X[i][j]);
        }
        printf("\nK = \n");
        for (j=0;j<ntot_states[i];j++){
            for (k=0;k<ntot_states[i];k++){
                printf("%f ",K[i][j][k]);
            }
            printf("\n");
        }
        printf("\nTrans = \n");
        for (j=0;j<ntot_states[i];j++){
	    printf("%d -> ",j);
            for (k=0;k<Num_trans[i][j];k++){
                printf("%d ",Trans[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    return;
}


double rand_exp(double k){
    double u;
    u = rand()/(RAND_MAX+1.0);
    return -log(1-u)/k;
}


void photophysics_t(double *t0, int *old_state, int *state0, double pos[3], double old_pos[3], double box[3], int **Nfs, int lamID, int **max_Nfs){
    double t=*t0,u;
    int i, move=0, new_state=*old_state, state=*old_state,last_fstate,comp,fluor,s_cnt;

    for(i=0;i<n_fluor_states[lamID];i++){
        Nfs[i][0]=0;
    }

    if(t>1){
        *t0=t-1;        //old_state is old_state and state0 is state0
        return;
    }

    for (i=0;i<3;i++){
        if(fabs(pos[i]-old_pos[i])>0.5*box[i]){
            move=1;
            break;
        }
    }
    if(move==1){
        new_state=0; // ground state S0
    }
    while(t<1){
        state=new_state;
        s_cnt=0;
        comp=0; // transition is complete
        fluor=0; // no emissions
        if (Num_trans[lamID][state]==0){ //no next state (for example photobleach)
            *t0=0;
            *old_state=state;
            *state0=state;
            return; //return Nfs as is.
        }

        for(i=0;i<Num_trans[lamID][state];i++){ //check infinite rate
            if (K[lamID][state][Trans[lamID][state][i]]<-0.5){
                comp=1;
                new_state=Trans[lamID][state][i];
                fluor=0;
                break;
            }
        }
        if(comp==0){//next state exists with a finite rate
            if (Num_trans[lamID][state]==1){//only one next state
                new_state=Trans[lamID][state][0];
            }
            else if (Num_trans[lamID][state]>1){
                u=rand()/(float)RAND_MAX;
                for (i=0;i<Num_trans[lamID][state];i++){
		    if (u<Prob[lamID][state][i]){
		        break;
		    }
                }
		new_state=Trans[lamID][state][i];
            }
            t+=rand_exp(K[lamID][state][new_state]);
            if(new_state>=n_states[lamID]){//emission
                Nfs[new_state-n_states[lamID]][0]++;
                last_fstate=new_state;
                new_state=Trans[lamID][new_state][0];
                fluor=1;
            }
        }
        if (t>1){
            if (fluor==1){
                Nfs[last_fstate-n_states[lamID]][0]--;
            }
            for (i=0;i<n_fluor_states[lamID];i++){
                if (Nfs[i][0]>max_Nfs[i][0]){
                    max_Nfs[i][0]=Nfs[i][0];
                }
            }
            *t0=t-1.0;
            *old_state=state;
            *state0=new_state;
            return; 
        }
    }
}

/* function: read_gro
 * ------------------
 * Reads a Gromacs structure file
 * 
 *   filename: Name of thr gromacs file
 *  ***fooIMG: Pointer containing image intensities
 *    ****psf: Pointer containing PSF.
 *   psf_type: 0 implies rotational symmetry and symmetry about object focal 
 *             plane 1 implies rotational symmetrt and not symmetric about 
 *             object focal plane
 */
 
void gro2spm(char filename[30], char outname[30]){
    FILE *f,*f2,*w;
    int cnt,cnt2, N, i, j,k, bb1=0, *cnt3, Nfluor=0;
    size_t width1,width2;
    char line[200], atom_name[6], line2[200], pos_str[9], *pend, *pend2, infile[70]; //Charcter start pointer
    char outfile[50];
    int **atom_idx, num_atoms[2],max_atoms[2], **state, **old_state;
    double ***pos, box[3], ***old_pos, **t0;

    cnt3=(int *)malloc(nlams*sizeof(int));
    atom_idx=(int **)malloc(nlams*sizeof(int *));
    for (i=0;i<2;i++){
        atom_idx[i]=(int *)malloc(50*sizeof(int));
        num_atoms[i]=0;
        max_atoms[i]=50;
    }

    ///////////////////////////////////////////////////////////////////////////
    /*                         READING INPUT FILE                            */
    ///////////////////////////////////////////////////////////////////////////
    f=fopen(filename,"r");
    if (f==NULL){
        printf("Input file '%s' does not exist\n",filename);
        return;
    }
    cnt=0;
    while ((fgets(line,200,f)) != NULL ){
        pend=strtok(line," \n");
        if(pend==NULL){
            continue;
        }
        if (cnt==0){ //first file
            f2 = fopen(line,"r");
            if (f2==NULL){
                printf("Input file %s does not exist\n",line);
                exit(-1);
            }
            cnt2=0;
            while (fgets(line2,200,f2)!=NULL){
                if (cnt2==1){
                    N=atoi(line2);
                }
                if (cnt2>1 && cnt2<N+2){
                    i=0;
                    while(line2[10+i]==' '){
                         i++; 
                    }
                    strncpy(atom_name,line2+10+i,5-i);
         
                    for (j=0;j<5-i;j++){
                        if (atom_name[j]==' '){
                            break;
                        }
                    }
                    atom_name[j]='\0';
                    for(i=0;i<nlams;i++){
                        for (j=0;j<nlam_atoms[i];j++){
                            if (strcmp(atom_name,lam_atoms[i][j])==0){//atomname matches
                                atom_idx[i][num_atoms[i]]=cnt2;
                                num_atoms[i]++;
                                if(num_atoms[i]==max_atoms[i]){
                                    atom_idx[i]=(int *)realloc(atom_idx[i],(max_atoms[i]+50)*sizeof(int));
                                    max_atoms[i]+=50;
                                }
                                break; //There is only one atom name so we should break.
                            }
                        }
                    }
                }
                else if (cnt2==N+2){

                    pend2=strtok(line2," \n");
                    box[0]=atof(pend2);
                    pend2=strtok(NULL," \n");
                    box[1]=atof(pend2);
                    pend2=strtok(NULL," \n");
                    box[2]=atof(pend2);
                }
                cnt2++;
            }
            fclose(f2);
            pos=(double ***)malloc(nlams*sizeof(double **));
            old_pos=(double ***)malloc(nlams*sizeof(double **));
            state=(int **)malloc(nlams*sizeof(int *));
            old_state=(int **)malloc(nlams*sizeof(int *));
            t0=(double **)malloc(nlams*sizeof(double *));
            
            for(i=0;i<nlams;i++){
                pos[i]=(double **)malloc(num_atoms[i]*sizeof(double *));
                old_pos[i]=(double **)malloc(num_atoms[i]*sizeof(double *));
                t0[i]=(double *)malloc(num_atoms[i]*sizeof(double));
                for(j=0;j<num_atoms[i];j++){
                    pos[i][j]=(double *)malloc(3*sizeof(double));
                    old_pos[i][j]=(double *)malloc(3*sizeof(double));
                    t0[i][j]=0.0;
                }
                state[i]=(int *)malloc(num_atoms[i]*sizeof(int));
                old_state[i]=(int *)malloc(num_atoms[i]*sizeof(int));
		for (j=0;j<num_atoms[i];j++){
		    state[i][j]=0;
		    old_state[i][j]=0;
		}
            }
        }
        ///read positions
        f2 = fopen(pend,"r");
        if (f2==NULL){
            printf("Input file %s does not exist\n",pend);
            exit(-1);
        }
        Nfluor=0;
        cnt2=0;
        for (i=0;i<nlams;i++){
            Nfluor+=num_atoms[i];
            cnt3[i]=0;
        }
        while (fgets(line2,200,f2)!=NULL){
            for (i=0;i<nlams;i++){
                if (atom_idx[i][cnt3[i]]==cnt2){
                    for(j=0;j<3;j++){
                        strncpy(pos_str,line2+20+j*8,8);
                        pos[i][cnt3[i]][j]=strtof(pos_str,&pend2);
                    }
                    cnt3[i]++;
                }
            }
            if (cnt2==N+2){
                pend2=strtok(line2," \n");
                box[0]=atof(pend2);
                pend2=strtok(NULL," \n");
                box[1]=atof(pend2);
                pend2=strtok(NULL," \n");
                box[2]=atof(pend2);
            }
 
            cnt2++;
        }
        if (cnt==0){//first file duplicate pos, state
            for (i=0;i<nlams;i++){
                for (j=0;j<num_atoms[i];j++){
                    for (k=0;k<3;k++){
                        old_pos[i][j][k]=pos[i][j][k];
                    }
                }
            }
        }

        i=0;j=0;k=0;
        while(strcmp2(line,".gro",k,k+4)){
            k++;
        }
        while(isdigit(line[k-i-1])){
           i++; 
        }
        line[k]='\0';
        sprintf(outfile,"%s%s.spm",outname,line+k-i);
        printf("Writing: %s         \r",outfile);
        w=fopen(outfile,"w");
        fprintf(w,"Nfluorophores %d\n",Nfluor);
        for (i=0;i<nlams;i++){
            for(j=0;j<num_atoms[i];j++){
                fprintf(w,"%10.3g%10.3g%10.3g%10s",pos[i][j][0],pos[i][j][1],pos[i][j][2],X[i][state[i][j]]);
                photophysics_t(&t0[i][j],&old_state[i][j],&state[i][j],pos[i][j],old_pos[i][j],box,Nfs[i],i,max_Nfs[i]);
                //old_state[i][j], state[i][j], Nfs are changed.
                fprintf(w,"%10s%10s",X[i][old_state[i][j]],X[i][state[i][j]]);
                fprintf(w,"%10.7f%5d",t0[i][j],n_fluor_states[i]);
                for (k=0;k<n_fluor_states[i];k++){
                    fprintf(w,"%10d%5d",Nfs[i][k][0],Nfs[i][k][1]);
                }
                fprintf(w,"\n");
                for (k=0;k<3;k++){
                    old_pos[i][j][k]=pos[i][j][k];
                }
            }
        }
        fprintf(w,"%f %f %f",box[0],box[1],box[2]);
        fclose(w);
	cnt++;
    }
    fclose(f);
    printf("\n");
    return;
}

void print_copyright(){
    printf("Copyright 2020,2021 SUBHAMOY MAHAJAN\n\n");
    printf("This file is part of InSilico Microscopy software\n\n");
    printf("InSilico Microsocpy is free software: you can redistribute it ");
    printf("and/or modify\nit under the terms of the GNU General Public ");
    printf("License as published by\nthe Free Software Foundation, either ");
    printf("version 3 of the License, or\n(at your option) any later ");
    printf("version.\n\nThis program is distributed in the hope that it ");
    printf("will be useful,\nbut WITHOUT ANY WARRANTY; without even the ");
    printf("implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR");
    printf(" PURPOSE.  See the\nGNU General Public License for more details.");
    printf("\n\nYou should have received a copy of the GNU General Public");
    printf(" License\nalong with this program.  If not, see <http://www.gn");
    printf("u.org/licenses/>.)\n");
}


int main(int argc, char* argv[]){
/*
   OPTIONS
   -f input filename (The gro file .gro)
   -o output filename (spm data file)
   -p parameters file (input)
   -pp photophysics file (input)
*///////////////////////////////Define Variables///////////////////////////////
    int i, j, k, l, val;
    char pp_file[50], infile[50], outfile[50], paramfile[50], filename[70], line[1200], *a,
         *varname, *pend, outmod[100];

/*  INTEGERS:
    i,j,k,l: are iteration variables
    val: generic value variable

    CHARACTERS:
    infile:  Input .gro structure file name.
    outfile:  Output image file name containing I(l',m'). {outname}_lam[0-9].dat.
    paramfile: Input parameters file.
    filename: Generic variable to store file names of different files.
    line: Generic variable to store lines read from a file.
    *a: Generic pointer to help retrieve doubles and int from a string
    *pend: Generic pointer to help retrieve doubles and int from a string

    FILE:
    *f: Generic file pointer (generated in code)
*//////////////////////////////////////////////////////////////////////////////

    for (i=0;i<10;i++){
        lams[i]=0;//initialize wavelengths to 0.
    }
    //read the arguments supplied
    for (i=0;i<argc;i++){
        if (strcmp(argv[i],"-silent")==0){ // Input gro file
            silent=1;
        }
    }
    if (silent==0) print_copyright();
    val=0; //number of arguments read
    for (i=2;i<argc;i+=2)
    {
        if (strcmp(argv[i-1],"-f")==0){ // Input gro file
            strcpy(infile,argv[i]);
            if (silent==0) printf("infile is %s\n",infile);
            val++;
        }
        else if (strcmp(argv[i-1],"-o")==0){ // Output image spm file
            strcpy(outfile,argv[i]);
            if (silent==0) printf("outfile is %s\n",outfile);
            val++;
        }
        else if (strcmp(argv[i-1],"-p")==0){ // Parameters file
            strcpy(paramfile,argv[i]);
            if (silent==0) printf("paramfile is %s\n",paramfile);
            val++;
        }
        else if (strcmp(argv[i-1],"-pp")==0){ // Photophysics file
            strcpy(pp_file,argv[i]);
            if (silent==0) printf("photophysic file is %s\n",pp_file);
            val++;
        }
    }
    if (val!=4){
        printf("Provide all arguments -f (input file), -o (output file), and");
        printf(" -p (parameters file -pp photophysics file)\n");
        return -1;
    }
    ///////////////////////////////////////////////////////////////////////////
    /*                         READING PARAMFILE                             */
    ///////////////////////////////////////////////////////////////////////////
    FILE *f;
    f=fopen(paramfile,"r");
    if (f==NULL){
        printf("Parameter file %s does not exist\n",paramfile);
        return -1;
    }
    while ((fgets(line,1200,f)) != NULL ){
        varname=strtok(line,"= ");
        if (varname[0]=='\0'){//empty line gives NULL variable names 
            continue;
        }
        a=strtok(NULL," =\n");

        if(strcmp2(varname,"lam",0,2)==0){

            if(strcmp2(varname,"_names",3,8)==0){
                val=varname[9]-'1';// Here val is saved as index of lambda
                nlam_atoms[val]=0;
                if (silent==0) printf("%s = ",varname);
                while(a!=NULL){
                    strcpy(lam_atoms[val][nlam_atoms[val]],a);
                    a=strtok(NULL," \n");
                    if (silent==0) printf("'%s' ",lam_atoms[val][nlam_atoms[val]]);
                    nlam_atoms[val]++;
                }
                if (silent==0) printf("\n");
                if (silent==0) printf("Number of atom names = %d\n",nlam_atoms[val]);
            }
            else if(strcmp2(varname,"_hue",3,6)==0){
                continue;
            }
            else if(strcmp2(varname,"_I0_",3,6)==0){
                continue;
            }
            else {
                val=varname[3]-'1';
                lams[val]=atoi(a);
                if (silent==0) printf("%s = %d \n",varname,lams[val]);
                if(val+1>nlams){
                    nlams=val+1;
                }
            }
        }
    }
    fclose(f);
    printf("Parameters File read\n");
    ///////////////////////////////////////////////////////////////////////////////

    srand((unsigned)time(NULL));
    ini_trans(pp_file);
    printf("Photophysics Initialized\n");
    gro2spm(infile,outfile);
    f=fopen(paramfile,"a");
    fprintf(f,"max_Nfs = ");
    for (i=0;i<nlams;i++){
        for (j=0;j<n_fluor_states[i];j++){
             fprintf(f,"%d ",max_Nfs[i][j][0]);
        }
        fprintf(f,"; ");
    }
    fprintf(f,"\n");
    fclose(f);
    return 0;
}
