/* Copyright 2020,2021 SUBHAMOY MAHAJAN

This program is part of in-silico-mocroscopy software package. 

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


/*
     opt_axis: Optical axis (n). 0, 1, and 2 for x, y, and z.
         nlam: Number of fluorophore types
   nlam_atoms: Number of atom names for each fluorophore type
    lam_atoms: Gromacs atom names for a fluorophore type. Stores upto 200
               atom names for upto 10 fluorophore types
           fs: is the multiplier multiplied to the MD structure coordinates or 
               1/fs is multiplied to lambda.
          lam: Wavelengths for each fluorophore type. 10 types supported.
         Nbox: is the number of small 3d boxes in each direction for the 
               coordinates.
         Npsf: The number of voxels in l', m', and n' direction for the PSF.
       MaxBox: Number of pixels in l and m directions for the image.
          pbc: Periodic boundary condition for each direction. 0 is false 1 is
               True
    lam_atoms: Contains all atoms names for all fluorophore types. 
    psfheader: Starting string of the file containing PSF data. 
           dx: dimensions of smallest voxels in each direction.
    focus_cor: Coordinate of focal plane.
       length: Dimensions of current gro simulation box in x, y, z.
         maxl: is the simulation (maximum w.r.t different time) box length in
               each direction; B*_x, B*_y, B*_z.
         Lpsf: length of PSF box length in each direction; P_l', P_m', P_n'
*/

int opt_axis, nlam, nlam_atoms[10], fs, lam[10], Npsf[3], Nbox[3], MaxBox[3],
    pbc[3]={0,0,0};
char lam_atoms[10][200][5], psfheader[100]; 
double dx[3], focus_cor, length[3], Lpsf[3], maxl[3], tsO=-1; 

/* function: min
 * -------------
 *  Finds miinimum of two doubles
 * 
 *  a: double number 1
 *  b: double number 2
 */
double min(double a, double b){
    if (a>b)
        return b;
    else
        return a;
}

/* function: read_psf
 * -----------------
   Reads the point spread function (PSF) data file. Each line should contain l',
   m', n' and PSF(l',m',n')

   filename: Filename containing PSF data
    ****psf: Pointer containing PSF of all fluorophore types.
     lam_id: ID of fluorophore type
   psf_type: 0 implies rotational symmetry and symmetry about object focal plane
             1 implies rotational symmetrt and not symmetric about object focal plane
 */

void read_psf(char filename[30],double  ****psf,int lam_id, int psf_type){
    FILE *f;
    char line[400], *pend;
    double x,y,z,I;
    int i,j,k,cnt=0;

    f=fopen(filename,"r");
    if (f==NULL){
        printf("PSF file %s does not exist\n",filename);
        exit(-1);
    }    
    while ((fgets(line,400,f)) != NULL ){
        if (line[0]=='#'){
            continue;
        }
        x=strtof(line,&pend); // l'
        y=strtof(pend,&pend); // m'
        z=strtof(pend,&pend); // n'
        I=strtof(pend,&pend); // PSF(l',m',n')
        
        // Ignoring lines were l', m', n' is not integral multiple 
        // of dl', dm', dn'. 
        if (x>0){
  	    i=(int)(x/dx[0]+1E-3);
        }
        else {
            i=(int)(x/dx[0]-1E-3);
        }
        if (cnt<10){
            if (fabs((double)i-x/dx[0])>dx[0]*0.1){
                printf("x: ignoring PSF line: %d not equal to %f. Difference is %f\n"
                       ,i,x/dx[0],fabs((double)i-x/dx[0]));
                cnt++;
                if (cnt==10){
                    printf("Similar output will not be shown as it repeated 10 times.\n");
                }
                continue;
            }
        }

        if (y>0){
            j=(int)(y/dx[1]+1E-3);
        }
        else{
            j=(int)(y/dx[1]-1E-3);
        }
        if (cnt<10){  
            if (fabs((double)j-y/dx[1])>dx[1]*0.1){
                printf("y: ignoring PSF line: %d not equal to %f. Difference is %f\n"
                       ,j,y/dx[1],fabs((double)j-y/dx[1]));
                cnt++;
                if (cnt==10){
                    printf("Similar output will not be shown as it repeated 10 times.\n");
                }
                continue;
            }
        }
        

        if (z>0){
            k=(int)(z/dx[2]+1E-3);
        }
        else{
            k=(int)(z/dx[2]-1E-3);
        }
        if (cnt<10){
            if (fabs((double)k-z/dx[2])>dx[2]*0.1){
                printf("z: ignoring PSF line: %d not equal to %f. Difference is %f\n"
                       ,k,z/dx[2],fabs((double)k-z/dx[2]));
                cnt++;
                if (cnt==10){
                    printf("Similar output will not be shown as it repeated 10 times.\n");
                }
                continue;
            }
        }

        if (psf_type==1){
            k+=(Npsf[2]-1)/2;
        }
            
        if (k>=Npsf[2]){
            break;
        }
        
        if ((i<Npsf[0])&&(j<Npsf[1])&&(k<Npsf[2])){
            psf[lam_id][i][j][k]=I;
            psf[lam_id][j][i][k]=I;
        }
    }
}

/* function: get_intensity
 * -----------------------
 *  Updates the image intensity I(l',m') contribution for a fluorophore particle
 *   located at (l_j,m_j,n_j) PSF(l'-l_j,m'-m_j,n'-n_j).
 *
 *    ****psf: Point spread function data 
 *     lam_id: ID of fluorophore type
 *  ***fooIMG: Image intensity I(l',m') data
 *        pos: Position of a fluorophore particle
 *   psf_type: 0 implies rotational symmetry and symmetry about object focal plane
 *             1 implies rotational symmetrt and not symmetric about object focal plane
 */

void get_intensity(double ****psf,int lam_id, double ***fooIMG, double pos[3], int psf_type){
// psf, fooIMG, Npsf, x, MaxBox are in image coordinates-- 
//      --->    no need to change direction based on opt_axis
// Nbox, pos, dx are in MD coordinates -- 
//      --->    need to change direction based on opt_axis
    int i, j, x[3], xid, yid, Mx, My;
    double dz;
    if (pos[(opt_axis+1)%3]>0){
        x[0]=(int)(pos[(opt_axis+1)%3]/dx[0]+0.5); // l axis
    }
    else {
        x[0]=(int)(pos[(opt_axis+1)%3]/dx[0]-0.5); // l axis
    }
    if (pos[(opt_axis+2)%3]>0){
        x[1]=(int)(pos[(opt_axis+2)%3]/dx[1]+0.5); // m axis
    }
    else {
        x[1]=(int)(pos[(opt_axis+2)%3]/dx[1]-0.5); // m axis
    }
   
    dz=focus_cor-pos[opt_axis];// n axis
    if (pbc[opt_axis]==1){//Check if PBC is applied in optical direction
        if (dz>length[opt_axis]*0.5){
            dz-=length[opt_axis]; // PBC considerations
        }
        else if (dz<-length[opt_axis]*0.5){
            dz+=length[opt_axis]; // PBC considerations
        }
    }
    if (psf_type==0){
        x[2]=(int)(fabs(dz)/dx[2]+0.5);
    }
    else if (psf_type==1){
        if (dz>0){
            x[2]=(int)(dz/dx[2]+0.5)+(Npsf[2]-1)/2;
        }
        else {
            x[2]=(int)(dz/dx[2]-0.5)+(Npsf[2]-1)/2;
        }
    }

    Mx=(MaxBox[0]-Nbox[0])/2;
    My=(MaxBox[1]-Nbox[1])/2;
    if ((x[2]<Npsf[2])&&(x[2]>=0)){
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
                if ((xid<Mx+Nbox[0])&&(yid<My+Nbox[1])&&(xid>=Mx)&&(yid>=My)){
                // To check if the image coordinate exists is in [0,Pl'] or 
                //  [0,Pm'].
                    fooIMG[lam_id][xid][yid]+= 
                                        psf[lam_id][abs(i)][abs(j)][x[2]];
                    //printf("xid=%d yid=%d, i=%d, j=%d, PSF=%f\n",xid,yid,i,j,psf[lam_id][abs(i)][abs(j)][x[2]]);
                }
            }
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
 *   psf_type: 0 implies rotational symmetry and symmetry about object focal plane
 *             1 implies rotational symmetrt and not symmetric about object focal plane
 */
void read_gro(char filename[30], double ***fooIMG, double ****psf, int psf_type){
    FILE *f;
    int cnt=0, N, i, j, bb1=0;
    char line[200], atom_name[6], pos_str[9], *pend; //Charcter start pointer
    double pos[3];
    f = fopen(filename,"r");
    if (f==NULL){
        printf("Input file %s does not exist\n",filename);
        exit(-1);
    }

    while (fgets(line,200,f)!=NULL){
        if (cnt==1){
            N=atoi(line);
        }
        else if (cnt>1 && cnt<N+2){
            for(i=0;i<3;i++){
                strncpy(pos_str,line+20+i*8,8);
                pos_str[8]='\0';
	        pos[i]=strtof(pos_str,&pend);
            }
            i=0;
            while(line[10+i]==' '){
                i++; }
            strncpy(atom_name,line+10+i,5-i);
            
            for (j=0;j<5-i;j++){
                if (atom_name[j]==' '){
                    break;
                }
            }
	    atom_name[j]='\0';

            for (i=0;i<nlam;i++){
                for (j=0;j<nlam_atoms[i];j++){
                    if (strcmp(atom_name,lam_atoms[i][j])==0){ //atomname matches
                        get_intensity(psf,i,fooIMG,pos,psf_type);//Update get_intensity
                        break; //There is only one atom name so we should break.
                    }
                }
            }
            /*For printing progress
            if (cnt%100==0){
                printf("%f%% done\n",(double)cnt*100.0/(double)N);
            }*/
        }
        cnt++;
    } 
    fclose(f);
    return;
}

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

int main(int argc, char* argv[] )
{
    print_copyright();
/*
   OPTIONS
   -f input filename (The gro file .gro)
   -o output filename (image data file)
   -p parameters file (input)
   -psf PSF header name (input)
*///////////////////////////////Define Variables////////////////////////////////
    int i, j, k, l, val, psf_type=0;
    char infile[50], outfile[50], paramfile[50], filename[70], line[1200], *a, 
         *varname, *pend, outmod[100];

/*  INTEGERS:
    i,j,k,l: are iteration variables
    val: generic value variable 
    psf_type: 0 implies rotational symmetry and symmetry about object focal plane (Gandy, default)
              1 implies rotational symmetry not symmetrical about object focal plane (Gibson and Lanni)
    currently both types are shift-invariant. For Depth-varying PSF different PSF should be read.

    CHARACTERS:
    infile: Input .gro structure file name. 
    outfile: Output image file name containing I(l',m'). {outname}_lam[0-9].dat.
    paramfile: Input parameters file. 
    filename: Generic variable to store file names of different files. 
    line: Generic variable to store lines read from a file. 
    *a: Generic pointer to help retrieve doubles and int from a string
    *pend: Generic pointer to help retrieve doubles and int from a string 
 
    FILE:
    *f: Generic file pointer (generated in code)
*//////////////////////////////////////////////////////////////////////////////

    for (i=0;i<10;i++){
        lam[i]=0;//initialize wavelengths to 0.
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
        else if (strcmp(argv[i-1],"-psf")==0){ // PSF header name file
            strcpy(psfheader,argv[i]);
            printf("psfheader is %s\n",psfheader);
            val++;
        }
    }
    if (val!=4){
        printf("Provide all arguments -f (input file), -o (output file), and");
        printf(" -p (parameters file) -psf (psf fileheader)\n");
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

        if(strcmp(varname,"fs")==0){
            fs=atoi(a);
            printf("%s = %d\n",varname,fs);
        }
        if(strcmp(varname,"tsO")==0){
            tsO=atof(a);
            printf("%s = %f\n",varname,tsO);
        }
        else if(strcmp(varname,"focus_cor")==0){
            focus_cor=atof(a);
            printf("%s = %f\n",varname,focus_cor);   
        }
        else if(strcmp(varname,"maxlen")==0){
            for (i=0;i<3;i++){
                maxl[i]=atof(a);
                a=strtok(NULL," \n");
            }
            printf("%s = [%f,%f,%f]\n",varname,maxl[0],maxl[1],maxl[2]);   
        }
        else if(strcmp(varname,"pbc")==0){
            while ((*a != '\0')&&(*a != '\n')){
                i=*a-'x';
                if ((i<=2)&&(i>=0)){
                   pbc[i]=1;
                }
                a++;
            }
            printf("%s = [%d,%d,%d] #1 is on, 0 is off\n", varname, pbc[0], 
                   pbc[1],pbc[2]);
        }
        else if(strcmp(varname,"opt_axis")==0){
            opt_axis=atoi(a);
            printf("%s = %d\n",varname,opt_axis);   
        }
        else if(strcmp2(varname,"lam",0,2)==0){
            if(strcmp2(varname,"_names",3,8)==0){
                val=varname[9]-'1';// Here val is saved as index of lambda
                nlam_atoms[val]=0;
                printf("%s = ",varname);
                while(a!=NULL){
                    strcpy(lam_atoms[val][nlam_atoms[val]],a);
                    a=strtok(NULL," \n");
                    printf("'%s' ",lam_atoms[val][nlam_atoms[val]]);
                    nlam_atoms[val]++;
                }
                printf("\n");
                printf("Number of atom names = %d\n",nlam_atoms[val]);
            }
            else if(strcmp2(varname,"_hue",3,6)==0){
                continue;
            }
            else if(strcmp2(varname,"_I0_",3,6)==0){
                continue;
            } 
            else {
                val=varname[3]-'1';
                lam[val]=atoi(a);   
                printf("%s = %d \n",varname,lam[val]);
                nlam++;
            }
        }
        else if(strcmp(varname,"dlmn")==0){
            for (i=0;i<3;i++){
                dx[i]=atof(a);
                a=strtok(NULL," \n");
            }
            if (fabs(dx[0]-dx[1])>1E-6){
                printf("Error! dx and dy should be the same.\n");
                return -1;
            }
            
            printf("%s = [%f,%f,%f]\n",varname,dx[0],dx[1],dx[2]);   
        }
        else if(strcmp(varname,"Plmn")==0){
            for (i=0;i<3;i++){
                Lpsf[i]=atof(a);
                a=strtok(NULL," \n");
            }
            printf("%s = [%f,%f,%f]\n",varname,Lpsf[0],Lpsf[1],Lpsf[2]);   
        }
        else if(strcmp(varname,"psf_type")==0){
            psf_type=atoi(a);
            printf("%s = %d\n",varname,psf_type);
        }
    }
    fclose(f); 
    ///////////////////////////////////////////////////////////////////////////
    // Read current box length; B_x, B_y, B_z in "length"
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
    ///////////////////////////////////////////////////////////////////////////

    for (i=0;i<3;i++){
        Npsf[i]=(int)(0.5*Lpsf[i]/dx[i]+1E-4)+1; 
        // Npsf: Number of PSF voxels in l', m', n'
        // length: is dimensions of the gro box in MD x, y, and z.
        // (opt_axis+i+1)%3 : is l, m, n directions, if i ix x, y, and z.
        // dx: voxel size dl, dm and dn
        // Number of voxels in MD simulation box in l, m, n
        Nbox[i]=(int)(length[(opt_axis+i+1)%3]/dx[i]+1E-3); 
    }
    if (psf_type==1){
        Npsf[2]=(int)(Lpsf[2]/dx[2]+1E-4)*2+1;
    }
    // Maximum number of voxels in largest MD simulation box in l and m. 
    MaxBox[0]=(int)(maxl[(opt_axis+1)%3]/dx[0]+1E-3); // Same as B*_l
    MaxBox[1]=(int)(maxl[(opt_axis+2)%3]/dx[1]+1E-3);// Same as B*_m

    printf("MaxBox: %d %d\n",MaxBox[0],MaxBox[1]);
    printf("Nbox: %d %d %d \n",Nbox[0],Nbox[1],Nbox[2]);
    printf("Npsf: %d %d %d \n", Npsf[0],Npsf[1],Npsf[2]);
    ///////////////////////////////////////////////////////////////////////////
    /*                          DEFINING IMAGE                               */
    ///////////////////////////////////////////////////////////////////////////
    double ***fooIMG=(double ***)malloc(nlam*sizeof(double **));
    for (i=0;i<nlam;i++)
    {
        fooIMG[i]=(double **)malloc(MaxBox[0]*sizeof(double *));
        for (j=0;j<MaxBox[0];j++)
        {
            fooIMG[i][j]=(double *)malloc(MaxBox[1]*sizeof(double));
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    /*                         INITIALIZING IMAGE                             */
    ////////////////////////////////////////////////////////////////////////////
    for (i=0;i<nlam;i++)
    {
        for (j=0;j<MaxBox[0];j++)
        {
            for (k=0;k<MaxBox[1];k++)
            {
                fooIMG[i][j][k]=0.0f;
            }
        }
        // Initialize white frame with -1.0. Since I(l',m') belongs to [0,1]
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
    ////////////////////////////////////////////////////////////////////////////
    /*                             DEFINING PSF                               */
    ////////////////////////////////////////////////////////////////////////////
    double ****psf=(double ****)malloc(nlam*sizeof(double ***));
    for (i=0;i<nlam;i++)
    {
        psf[i]=(double ***)malloc(Npsf[0]*sizeof(double **));
        for (j=0;j<Npsf[0];j++)
        {
            psf[i][j]=(double **)malloc(Npsf[1]*sizeof(double *));
            for (k=0;k<Npsf[1];k++)
            {
                psf[i][j][k]=(double *)malloc(Npsf[2]*sizeof(double));
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    /*                            INITIALIZING PSF                           */
    ///////////////////////////////////////////////////////////////////////////
    for (i=0;i<nlam;i++){
        for (j=0;j<Npsf[0];j++){
            for (k=0;k<Npsf[1];k++){
                for (l=0;l<Npsf[2];l++){
                    psf[i][j][k][l]=0.0f;
                }
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    printf("Initialization of PSF Done\n");   
    //Read PSF
    for(i=0;i<nlam;i++){
        if (psf_type==0){
            sprintf(filename,"%s_lam%d_fs%d.dat",psfheader,lam[i],fs);
        }
        else if (psf_type==1) {
            sprintf(filename,"%s_tsO%g_lam%d_fs%d.dat",psfheader,tsO,lam[i],fs);
        }
        read_psf(filename,psf,i,psf_type);
        printf("%s: PSF for %d nm read.\n",filename,lam[i]);
    }
    printf("PSF read\nReading Gro\n");

    // Read gro and calculate I(l',m')
    read_gro(infile,fooIMG,psf,psf_type);
    printf("Analyzing Gro done\n");

    // Writing I(l',m')
    printf("Writing data files\n");
    for (l=0;l<nlam;l++)
    {
        FILE *f;
        if (psf_type==0){
            sprintf(outmod,"%s_lam%d_fs%d.dat",outfile,lam[l],fs);
        }
        else if (psf_type==1){
            sprintf(outmod,"%s_tsO%g_lam%d_fs%d.dat",outfile,tsO,lam[l],fs);
        }
        printf("output file is %s\n",outmod);
        f=fopen(outmod,"w");
        fprintf(f,"#Lam= %d dx= %f,%f,%f fs= %d\n",lam[l],dx[0],dx[1],dx[2],fs);

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
