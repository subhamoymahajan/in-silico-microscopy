/* Copyright 2020 SUBHAMOY MAHAJAN

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.) */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#define PI 3.141592653589793

int opt_axis,nlam,nlam_atoms[10],fs,lam[10],Npsf[3],Nbox[3],MaxBox[3],pbc[3]={0,0,0};
char lam_atoms[10][200][5],psfheader[100]; // psftype removed.
float  dx[3],focus_cor,length[3],Lpsf[3],NA,maxl[3]; 
float min(float a, float b){
    if (a>b)
        return b;
    else
        return a;
}

void read_psf(char filename[30],float  *****psf,int lam_id){
    FILE *f;
    char line[200];
    char *pend;
    char *pstart;
    float x,y,z,I;
    int i,j,k;//,cnt=0;
    f=fopen(filename,"r");
    if (f==NULL){
        printf("PSF file %s does not exist\n",filename);
        exit(-1);
    }    
    while ((fgets(line,200,f)) != NULL ){
        if (line[0]=='#'){
           continue;
        }
        x=strtof(line,&pend);
        y=strtof(pend,&pend);
        z=strtof(pend,&pend);
        I=strtof(pend,&pend);
        i=(int)(x/dx[0]+1E-3);
        if (abs((float)i-x/dx[0])>1E-3){
            printf("x: ignoring PSF line: %d not equal to %f. Difference is %f\n",i,x/dx[0],abs((float)i-x/dx[0]));
            continue;
        }
        j=(int)(y/dx[1]+1E-3);
        if (abs((float)j-y/dx[1])>1E-3){
            printf("y: ignoring PSF line: %d not equal to %f. Difference is %f\n",j,y/dx[1],abs((float)j-y/dx[1]));
            continue;
        }
        k=(int)(z/dx[2]+1E-3);
        if (abs((float)k-z/dx[2])>1E-3){
            printf("z: ignoring PSF line: %d not equal to %f. Difference is %f\n",k,z/dx[2],abs((float)k-z/dx[2]));
            continue;
        }
            
        if (k>=Npsf[2]){
            break;
        }
        if ((i<Npsf[0])&&(j<Npsf[1])&&(k<Npsf[2])){
            (*psf)[lam_id][i][j][k]=I;
            (*psf)[lam_id][j][i][k]=I;
        }
    }
}

void get_intensity(float *****psf,int lam_id, float ****fooIMG, float pos[3]){
//Caution: Make sure dx is same for python code and C code.
// psf, fooIMG, Npsf, x, MaxBox are in image coordinates-- no need to change direction based on opt_axis
// Nbox, pos, dx are in MD coordinates -- need to change direction based on opt_axis
    int i,j,a,x[3],xid,yid,Mx,My;
    float tempxy,dz;
    x[0]=(int)(pos[(opt_axis+1)%3]/dx[0]+1E-3); // l axis
    x[1]=(int)(pos[(opt_axis+2)%3]/dx[1]+1E-3); // m axis
    dz=fabs(pos[opt_axis]-focus_cor);// n axis
    tempxy=pos[(opt_axis+1)%3]/dx[0]-(float)x[0]; //To ensure atom is moved to the closed grid point 
    if (tempxy>0.5){
         x[0]++;
    }
    tempxy=pos[(opt_axis+2)%3]/dx[1]-(float)x[1]; //To ensure atom is moved to the closed grid point.
    if (tempxy>0.5){
         x[1]++;
    }
    if (pbc[opt_axis]==1){//Check if PBC is applied in optical direction
        if (dz>length[opt_axis]*0.5){
             dz=fabs(dz-length[opt_axis]); // PBC considerations
        }
    }
    x[2]=(int)(dz/dx[2]+1E-3);
    Mx=(MaxBox[0]-Nbox[0])/2;
    My=(MaxBox[1]-Nbox[1])/2;

    if (abs(x[2])*2<Npsf[2]){
        for (i=-Npsf[0]+1;i<Npsf[0];i++){
            for (j=-Npsf[1]+1;j<Npsf[1];j++){
                xid=i+x[0];
                yid=j+x[1];
                //PBC considerations
                if (pbc[(opt_axis+1)%3]==1){
                   xid=(xid+2*Nbox[0])%Nbox[0];
                }
                if (pbc[(opt_axis+2)%3]==1){
                   yid=(yid+2*Nbox[1])%Nbox[1];
                }
                // Centering the image in white background
                xid+=Mx;
                yid+=My;
                if ((xid<Mx+Nbox[0])&&(yid<My+Nbox[1])&&(xid>=Mx)&&(yid>=My)){// To check if the image coordinate exists when PBC is not applied.
                   (*fooIMG)[lam_id][xid][yid]+=(*psf)[lam_id][abs(i)][abs(j)][abs(x[2])];
                }
            }
        }
    }
}
void get_atom_name(char line[70], char atom_name[6]){
    int i,cnt=0;
    for(i=10;i<15;i++){
        if(line[i]!=' '){
            atom_name[cnt]=line[i];
            cnt++;
        }
    }
    atom_name[cnt]='\0';
}


void read_gro(char filename[30], float ****fooIMG, float *****psf){
    FILE *f;
    int cnt=0, N=10,i,j,acnt;
    char line[70], atom_name[6], atomcor[9], *pend, *pstart; //Charcter start pointer
    float pos[3];
    f = fopen(filename,"r");
    if (f==NULL){
       printf("Input file %s does not exist\n",filename);
       exit(-1);
    }
    char eightspace[9]="        ";
    while (cnt<N+3){
        fgets(line,70,f);
        if (cnt==1){
            N=atoi(line);
        }
        else if (cnt>1 && cnt<N+2){
            /*pstart=(char*)&line+20;
            pos[0]=strtof(pstart,&pend);
            pos[1]=strtof(pend,&pend);
            pos[2]=strtof(pend,&pend);*/
            acnt=20;
            for (i=0;i<3;i++){
                strcpy(atomcor,eightspace); //initialize to 8 spaces
                for (j=0;j<8;j++){
                    if((line[acnt]!='\n')&&(line[acnt]!='\0')){
                    	atomcor[j]=line[acnt];    
                    	acnt++;
                    }
                } 
                atomcor[j]='\0';
                pos[i]=strtof(atomcor,&pend);
            }//Verify this. new addition
            get_atom_name(line,atom_name);
            for (i=0;i<nlam;i++){
                for (j=0;j<nlam_atoms[i];j++){
                   if (strcmp(atom_name,lam_atoms[i][j])==0){ // atomname matches 
                       get_intensity(psf,i,fooIMG,pos);//Update get_intensity
                       break; // newly added. There is only one atom name so we can break.
                   }
                }
            }
/*For printing progress
            if (cnt%100==0){
                printf("%f%% done\n",(float)cnt*100.0/(float)N);
            }*/
        }
        cnt+=1;
    } 
    fclose(f);
    return;
}

char* find_var_name(char line[1200], char delim, char varname[10]){
// finds variable name from a provided line read from the parameters file.
    int i,cnt=0;
    char *a;
    varname[0]='\0'; 
    for(i=0;i<1200;i++){
        if((line[i]==' ')&&(cnt==0)){//ignore starting spaces
           continue;
        }
        else if(line[i]=='\n'){//line ends: stop reading
           break;
        }
        else if(line[i]=='='){//variable name exists only before =
            varname[cnt]='\0';
            break;
        }
        else if(line[i]!=' '){//add characters to variable name.
           varname[cnt]=line[i];
           cnt+=1;
        }
    }
     
    a=&(line[i+1]); // returns a pointer which points after '=' symbol.
    return a;
}
int strcmp2(char *a, char *b, int n1,int n2){
//Compares string a[n1:n2] with string b[0:n2-n1] 
    int i, res=0; //res=0 is success
    
    for(i=n1;i<=n2;i++){
       res=abs(a[i] - b[i-n1]) ;
    }
    return res;
}

int get_lamnames( char *a, int val, char lam_atoms[10][200][5]){
//stores atom names for each lambda and returns number of atom names for the lambda.

    char element[6]; //each atom name in gro file is 5 characters long
    int ele_idx=0, atom_idx=0,i=0; 
    element[0]='\0';
  
    while(*(a+i) != '\0' ){
        if((*(a+i) ==' ')&&(ele_idx==0)){ // ignore strating spaces
            i=i+1;
        }
        else if ((*(a+i)==' ')||(*(a+i)=='\n')){ // the atom name has ended
            element[ele_idx]='\0';
            ele_idx=0; //reinitialize for new atom name
            strcpy(lam_atoms[val][atom_idx],element);
            atom_idx+=1;
            i=i+1;
        }
        else {
            element[ele_idx]=*(a+i);
            ele_idx+=1;
            i=i+1;
        }
          
    }
    return atom_idx;
}
void print_copyright(){
    printf("Copyright 2020 SUBHAMOY MAHAJAN\n\nThis file is part of InSilico Microscopy software\n\nInSilico Microsocpy is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License\nalong with this program.  If not, see <http://www.gnu.org/licenses/>.)\n");
} 

int main(int argc, char* argv[] )
{
    print_copyright();
/*
   OPTIONS
   -f input filename (The gro file .gro)
   -o output filename (image data file)
   -p parameters file (input)
*/
/////////////////////////////////////Define Variables////////////////////////////////////////////
    int i,j,k,l,val;
    char infile[50], outfile[50],paramfile[50],filename[70],line[1200],*a,varname[10],*pend,outmod[100];
    float maxsize[2]={0.0,0.0};
/*
    INTEGERS:
    i,j,k: are iteration variables
    fs: is the multiplier multiplied to the MD structure coordinates of 1/fs is multiplied to lambda (for computational efficiency). (Default value 800)
    Nbox: is the number of small 3d boxes in each direction for the coordinates. (generated in code)
    Npsf: is the number of small boxes in each direction for the PSF file. (generated in code)
    lam: stores upto 10 different lambda values and in tested only for 3 lambdas. (Values must be present in Parameter file)
    opt_axis: 0 for x axis, 1 for y axis 2 for z axis. (Default value 2)
    timesteps: Total number of timesteps in the simulation
    val: generic value variable 

    FLOATS:
    length: is the simulation (maximum w.r.t different time) box length in each direction (Values must be present in the Parameter file)
    focus_cor: is the coordinate of focal plane (default value is length/2.0)
    dx: dimensions of smallest voxels in each direction.
    Lpsf: length of PSF box length in each direction (Values must be present in the parmeter file)

    CHARACTERS:
    infile: Input .gro structure file name (Must be GRO). Reading gro file not extensively tested. (Must be provided as an option)
    outfile: Output image file name. {outname}_lam[0-9].dat files will be created depending on number of lambdas (Must be provided as an option)
    paramfile: Input parameters file. It must be provided as some parameters do not have default values. (Must be provided as an option)
    boxfile: Input file containing box dimensions for each time.
    filename: Generic variable to store file names of different files. (generated in code)
    line: Generic variable to store lines read from a file. (generated in code)
    *a: Generic pointer to help retrieve floats and int from a string (generated in code)
    *pend: Generic pointer to help retrieve floats and int from a string (generated in code)
    lam_atoms: the atom names associated with each lambda. stores upto 200 different atoms names for upto 10 different lambdas.
 
    FILE:
    *f: Generic file pointer (generated in code)
*/
/////////////////////////////////////////////////////////////////////////////////////////////////
    for (i=0;i<10;i++){
        lam[i]=0;//initialize wavelengths to 0.
        for (j=0;j<200;j++){
            lam_atoms[i][j][0]='\0';//initialize atom names to null strings.
        }
    }
    //read the arguments supplied
    val=0; //number of arguments read
    for (i=2;i<argc;i+=2)
    {
       if (strcmp(argv[i-1],"-f")==0){ // Input gro file
           strcpy(infile,argv[i]);
           printf("infile is %s\n",infile);
           val++;
       }
       else if (strcmp(argv[i-1],"-o")==0){ // Output image data file
           strcpy(outfile,argv[i]);
           printf("outfile is %s\n",outfile);
           val++;
       }
       else if (strcmp(argv[i-1],"-p")==0){ // Parameters file
           strcpy(paramfile,argv[i]);
           printf("paramfile is %s\n",paramfile);
           val++;
       }
    }
    if (val!=3){
       printf("Provide all arguments -f (input file), -o (output file), and -p (parameters file)\n");
       return -1;
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    /*                                    READING PARAMFILE                                   */
    //////////////////////////////////////////////////////////////////////////////////////////
    FILE *f;
    f=fopen(paramfile,"r");
    if (f==NULL){
       printf("Parameter file %s does not exist\n",paramfile);
       return -1;
    }
    while ((fgets(line,1200,f)) != NULL ){
       a=find_var_name(line,'=',varname);
       if (varname[0]=='\0'){//empty line gives NULL variable names 
           continue;
       }
       printf("variable '%s' = ",varname);

       if(strcmp(varname,"f")==0){
           fs=atoi(a);
           printf("%d\n",fs);
       }
       else if(strcmp(varname,"focus_cor")==0){
           focus_cor=strtof(a,&pend);
           printf("%f\n",focus_cor);   
       }
       else if(strcmp(varname,"NA")==0){
           NA=strtof(a,&pend);
           printf("%f\n",NA);
       }
       else if(strcmp(varname,"maxlen")==0){
           maxl[0]=strtof(a,&pend);
           maxl[1]=strtof(pend,&pend);
           maxl[2]=strtof(pend,&pend);
           printf("[%f,%f,%f]\n",maxl[0],maxl[1],maxl[2]);   
       }
       else if(strcmp(varname,"pbc")==0){
           
           while ((a != NULL)&&(*a != '\n')){
              if (*a == 'x'){
                 pbc[0]=1;
              }
              else if (*a == 'y'){
                 pbc[1]=1;
              }
              else if (*a == 'z'){
                 pbc[2]=1;
              }
              a++;
           }
           printf("[%d,%d,%d] #1 is on, 0 is off\n",pbc[0],pbc[1],pbc[2]);
       }
       else if(strcmp(varname,"psfheader")==0){
           i=0;
           j=0;
           while((a[i]!='\0')&&(a[i]!='\n')){
               if (a[i]==' '){
                   if (j>0){
                       break; //psfreader cannot have space.
                   }
               }
               else {
                   psfheader[j]=a[i];
                   j++;
               }
               i++;
           }
           psfheader[j]='\0';
           printf("%s\n",psfheader);
       }
       else if(strcmp(varname,"opt_axis")==0){
           opt_axis=atoi(a);
           printf("%d\n",opt_axis);   
       }
       else if(strcmp2(varname,"lam",0,2)==0){
           if(strcmp2(varname,"_names",3,8)==0){
               val=varname[9]-'0'-1;// Here val is saved as index of lambda
               nlam_atoms[val]=get_lamnames(a,val,lam_atoms);
               for (i=0;i<200;i++){
                  if(lam_atoms[val][i][0]!='\0'){ 
                      printf("%s ",lam_atoms[val][i]);
                  }
               }
               printf("\n");
           }
           else {
               val=varname[3]-'0'-1;
               lam[val]=atoi(a);   
               printf("%d \n",lam[val]);
               nlam++;
           }
       }
       else if(strcmp(varname,"dx")==0){
           dx[0]=strtof(a,&pend);
           dx[1]=strtof(pend,&pend);
           dx[2]=strtof(pend,&pend);
           if (fabs(dx[0]-dx[1])>1E-6){
               printf("Error! dx and dy should be the same.\n");
               return -1;
           }
           
           printf("[%f,%f,%f]\n",dx[0],dx[1],dx[2]);   
       }
       else if(strcmp(varname,"Lpsf")==0){
           Lpsf[0]=strtof(a,&pend);
           Lpsf[1]=strtof(pend,&pend);
           Lpsf[2]=strtof(pend,&pend);
           printf("[%f,%f,%f]\n",Lpsf[0],Lpsf[1],Lpsf[2]);   
       }
       else
           printf("\n");
    }
    fclose(f); 
    //////////////////////////////////////////////////////////////////////////////////////////
    //reading current box length
    f=fopen(infile,"r");
    if (f==NULL){
       printf("Input file %s does not exist\n",infile);
       return -1;
    }
    i=0;
    j=0;
    while ((fgets(line,1200,f)) != NULL ){
        if (i==1){
            j=atoi(line); //Number of atoms
        }
        else if (i==2+j){
            length[0]=strtof(line,&pend);
            length[1]=strtof(pend,&pend);
            length[2]=strtof(pend,&pend);
            printf("Box length is [%f,%f,%f]\n",length[0],length[1],length[2]);
        }
        i++; 
    }
    fclose(f);    
    ////////

     
    for (i=0;i<3;i++){
        Npsf[i]=(int)(0.5*Lpsf[i]/dx[i]+1E-4)+1;
//length is dimensions of the gro box in MD x, y, and z (which can be different for image x and y depending on the opt_axis)
//dx [image x , image y, perpendicular to the image ] calculated for PSF.
        Nbox[i]=(int)(length[(opt_axis+i+1)%3]/dx[i]+1E-3); // Nbox [ image x, image y, and perpendicular to image]
    }
    //Max dimensions of an image is the max dimensions of the MD box in directions other than the focus coordinate
    MaxBox[0]=(int)(maxl[(opt_axis+1)%3]/dx[0]+1E-3); // Max image x
    MaxBox[1]=(int)(maxl[(opt_axis+2)%3]/dx[1]+1E-3);// Max image y

    printf("MaxBox: %d %d\n",MaxBox[0],MaxBox[1]);
    printf("Nbox: %d %d %d \n",Nbox[0],Nbox[1],Nbox[2]);
    printf("Npsf: %d %d %d \n", Npsf[0],Npsf[1],Npsf[2]);
    //////////////////////////////////////////////////////////////////////////////////////////
    /*                                    DEFINING IMAGE                                    */
    //////////////////////////////////////////////////////////////////////////////////////////
    float ***fooIMG=(float ***)malloc(nlam*sizeof(float **));
    for (i=0;i<nlam;i++)
    {
        fooIMG[i]=(float **)malloc(MaxBox[0]*sizeof(float *));
        for (j=0;j<MaxBox[0];j++)
        {
            fooIMG[i][j]=(float *)malloc(MaxBox[1]*sizeof(float));
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    /*                                 INITIALIZING IMAGE                                   */
    //////////////////////////////////////////////////////////////////////////////////////////
    for (i=0;i<nlam;i++)
    {
        for (j=0;j<MaxBox[0];j++)
        {
            for (k=0;k<MaxBox[1];k++)
            {
                fooIMG[i][j][k]=0.0f;
            }
        }
        //Newly added
        for (j=0;j<MaxBox[0];j++)
        {
            for (k=0;k<(MaxBox[1]-Nbox[1])/2;k++)
            {
                fooIMG[i][j][k]=-1.0;
            }
            for (k=(MaxBox[1]-Nbox[1])/2+Nbox[1];k<MaxBox[1];k++)
            {
                fooIMG[i][j][k]=-1.0;
            }
        }
        for (k=0;k<MaxBox[1];k++)
        {
            for (j=0;j<(MaxBox[0]-Nbox[0])/2;j++)
            {
                fooIMG[i][j][k]=-1.0;
            }
            for (j=(MaxBox[0]-Nbox[0])/2+Nbox[0];j<MaxBox[0];j++)
            {
                fooIMG[i][j][k]=-1.0;
            }
        }
    }
    printf("Initialization of Image Done\n");    
    //////////////////////////////////////////////////////////////////////////////////////////
    /*                                    DEFINING PSF                                      */
    //////////////////////////////////////////////////////////////////////////////////////////
    float ****psf=(float ****)malloc(nlam*sizeof(float ***));
    for (i=0;i<nlam;i++)
    {
        psf[i]=(float ***)malloc(Npsf[0]*sizeof(float **));
        for (j=0;j<Npsf[0];j++)
        {
            psf[i][j]=(float **)malloc(Npsf[1]*sizeof(float *));
            for (k=0;k<Npsf[1];k++)
            {
                 psf[i][j][k]=(float *)malloc(Npsf[2]*sizeof(float));
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    /*                                INITIALIZING PSF                                      */
    //////////////////////////////////////////////////////////////////////////////////////////
   for (i=0;i<nlam;i++){
       for (j=0;j<Npsf[0];j++){
           for (k=0;k<Npsf[1];k++){
               for (l=0;l<Npsf[2];l++){
                   psf[i][j][k][l]=0.0f;
               }
           }
       }
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////////
   printf("Initialization of PSF Done\n");   
   for(i=0;i<nlam;i++){
        sprintf(filename,"%s_lam%d_fs%d.dat",psfheader,lam[i],fs);
        read_psf(filename,&psf,i);
        printf("%s: PSF for %d nm read.\n",filename,lam[i]);
   }
   printf("PSF read\nReading Gro\n");
   //Centering the Image
   read_gro(infile,&fooIMG,&psf);
   printf("Analyzing Gro done\n");
   printf("Writing data files\n");
   for (l=0;l<nlam;l++)
   {
       FILE *f;
       sprintf(outmod,"%s_lam%d_fs%d.dat",outfile,lam[l],fs);
       printf("output file is %s\n",outmod);
       f=fopen(outmod,"w");
       fprintf(f,"#NA= %f Lam= %d dx= %f,%f,%f fs= %d\n",NA,lam[l],dx[0],dx[1],dx[2],fs);

       for (i=0;i<MaxBox[0];i++)
       {
           for (j=0;j<MaxBox[1];j++)
           {
               fprintf(f,"%f ",fooIMG[l][i][j]);
           }
           fprintf(f,"\n");
       }
       fclose(f);
       f=NULL;
    }
 return 0;
}
