#include "b1_timing.hpp"

int main(int nArgs, char* Args[])
{
    const int Infinity=-1;
    using NumT = bertini::dbl;

    std::chrono::time_point<std::chrono::high_resolution_clock> ct1_b1, ct2_b1, ct1_b2, ct2_b2, for_start, for_end;//initiate timing parameters
    std::chrono::duration<double> elapsed_seconds_b1, elapsed_seconds_b2, for_laps;

    FILE *fpA_b1, *fpx_b1, *fpb_b1, *fpelem_b1;//initiate file pointers for b1
    // std::ofstream fpA_b2, fpx_b2, fpb_b2;
	if(nArgs<4)
	{
		int prec_b1=1024, size=10, iter=1;//test parameters: prec_b1= bit precision, size=square matrix size, iter=number of repetition
		std::cout << "Not enough argument entered.\n Program will be executed with default parameters.\n bit_precision="<<prec_b1<<", size="<<size<<"iterations="<<iter<<" \n";
	}
	else	
	{
		int prec_b1=atoi(Args[1]), size=atoi(Args[2]), iter=atoi(Args[3]);
	}
    initMP(prec_b1);

    int prec_b2=log10(prec_b1)/log2(prec_b1)*prec_b1;//convert precision in bits for b1 in precision in digits for b2
    // std::cout<<prec_b2<< "s\n";//uncomment to check number of digits

    double tol=5e-324, largeChange=1e16;//parameters for some solvers in b1
    mpf_t mpftol,  mpflargeChange;

    mpf_init2(mpftol,prec_b1);//different initialization of parameters for some solvers in b1 /*needs to be checked*/
    mpf_init2(mpflargeChange,prec_b1);
    mpf_set_d(mpftol,tol);
    mpf_set_d(mpflargeChange,largeChange);
    double cond_num,  norm_A,  norm_A_inv;//auxiliary variables for some solvers in b1

    uint usize=(int)size;//convert size to unsigned integer /*needed for RandomOfUnits*/

    bertini::mpfr_float::default_precision(prec_b2);//set mp precision for b2

    bertini::Vec<bertini::complex> x_b2;//initialize solution vector in b2

    mat_mp A_b1;//initialize and set precision for b1 matrix and vectors
    vec_mp x_b1, b_b1;
    init_mat_mp2(A_b1,size,size,prec_b1);
    init_vec_mp2(x_b1,size,prec_b1);
    init_vec_mp2(b_b1,size,prec_b1);

    bool DontAlignCols=false;//IO parameters for Eigen IOFormat
    int sprec=prec_b2;

    Eigen::IOFormat CommaInitFmt(sprec, DontAlignCols, ", ", ", ", "", "", " << ", ";");//different Eigen IOFormats
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    Eigen::IOFormat OctaveFmt(sprec, 0, ", ", ";\n", "", "", "[", "]");
    Eigen::IOFormat MatlabFmtInline(sprec, 0, ", ", ";", "", "", "[", "]");//Matlab format one line
    Eigen::IOFormat HeavyFmt(sprec, 0, ", ", ";\n", "[", "]", "[", "]");
    std::string sep = "\n----------------------------------------\n";

    for_start=std::chrono::high_resolution_clock::now();//start count for overall time

    int tmp=0;//for-loop varaibles
    int i,j,k;

    // bertini::Mat<bertini::complex> A_b2=bertini::RandomOfUnits<bertini::complex>(usize,usize).eval();//define random matrix and vectors in b2 /*Note values are in the unit circle but the Matrix is NOT orthonormal*/
    // bertini::Vec<bertini::complex> b_b2=bertini::RandomOfUnits<bertini::complex>(usize).eval();
    bertini::Mat<bertini::complex> A_b2(size,size);//initialize matrix and vector in b2
    bertini::Vec<bertini::complex> b_b2(size);


	for (i=0;i<iter;i++)//start iterations
		{
		make_matrix_random_mp(A_b1,size, size, prec_b1);//define random matrix and vectors in b1 /*Note values are in the unit circle and the Matrix is orthonormal*/
		make_vec_random_mp(b_b1,size);
    /*testing*/
    // _point_mp elem_b1=b_b1[0];
    // _point_mp elem_b1=1;
    // _point_mp elem_b1=b_b1[0]->r;
    // comp_mp elem;
    // init_mp2(elem, prec_b1)
    // set_mp(elem,&b_b1->coord[0]);
    // fpelem_b1 = fopen ("elemb_b1.txt" , "a+");
		// print_Matlab_mp(fpelem_b1, sprec, elem);
    // fclose(fpelem_b1);
    // clear_mp(elem);

    // bertini::Mat<bertini::complex> A_b2=bertini::RandomOfUnits<bertini::complex>(usize,usize).eval();//define random matrix and vectors in b2 /*Note values are in the unit circle but the Matrix is NOT orthonormal*/
    // bertini::Vec<bertini::complex> b_b2=bertini::RandomOfUnits<bertini::complex>(usize).eval();
    /*testing*/
    // bertini::complex elem_b2=b_b2[0];
    // bertini::complex elem_b2=1;
    // b_b2[0]=1;
    // std::cin>>b_b2[0];
    // std::cin>>b_b1[0];
    // elem_b2=b_b1[0];

/*the following routine takes the vector b_b1 and the Matrix A_b1 defined in b1 and uses their elements to generate the corresponding vector and matrix in b2 format b_b2 and A_b2*/
    comp_mp vec_elem, mat_elem;
    init_mp2(vec_elem, prec_b1);
    init_mp2(mat_elem, prec_b1);
    char *ch1 = NULL, *ch2 = NULL, buffer [2*prec_b2+10];
    long e1,e2;
    for(j=0;j<size;j++)
      {
        for(k=0;k<size;k++)
          {
            set_mp(mat_elem,&A_b1->entry[j][k]);
            if (mpfr_number_p(mat_elem->r) && mpfr_number_p(mat_elem->i))
            {
              ch1 = mpf_get_str(NULL, &e1, 10, prec_b2, mat_elem->r);
              ch2 = mpf_get_str(NULL, &e2, 10, prec_b2, mat_elem->i);

              if (ch1[0] != '-')
                if (ch2[0] != '-')
                  sprintf(buffer, "(0.%se%ld,0.%se%ld)", ch1, e1, ch2, e2);
                else
                  sprintf(buffer, "(0.%se%ld,-0.%se%ld)", ch1, e1, &ch2[1], e2);
              else
                if (ch2[0] != '-')
                  sprintf(buffer, "(-0.%se%ld,0.%se%ld)", &ch1[1], e1, ch2, e2);
                else
                  sprintf(buffer, "(-0.%se%ld,-0.%se%ld)", &ch1[1], e1, &ch2[1], e2);

              free(ch1);
              free(ch2);
            }
            else
              sprintf(buffer, "(NaN,NaN)");

            std::istringstream buffer_iss (buffer);
            buffer_iss>>A_b2(j,k);

            if(k==0)
            {
              set_mp(vec_elem,&b_b1->coord[j]);
              if (mpfr_number_p(vec_elem->r) && mpfr_number_p(vec_elem->i))
              {
                ch1 = mpf_get_str(NULL, &e1, 10, prec_b2, vec_elem->r);
                ch2 = mpf_get_str(NULL, &e2, 10, prec_b2, vec_elem->i);

                if (ch1[0] != '-')
                  if (ch2[0] != '-')
                    sprintf(buffer, "(0.%se%ld,0.%se%ld)", ch1, e1, ch2, e2);
                  else
                    sprintf(buffer, "(0.%se%ld,-0.%se%ld)", ch1, e1, &ch2[1], e2);
                else
                  if (ch2[0] != '-')
                    sprintf(buffer, "(-0.%se%ld,0.%se%ld)", &ch1[1], e1, ch2, e2);
                  else
                    sprintf(buffer, "(-0.%se%ld,-0.%se%ld)", &ch1[1], e1, &ch2[1], e2);

                free(ch1);
                free(ch2);
              }
              else
                sprintf(buffer, "(NaN,NaN)");

              // std::cout<<buffer;
              std::istringstream buffer_iss (buffer);
              buffer_iss>>b_b2[j];
            }

          }
      }
    clear_mp(vec_elem);
    clear_mp(mat_elem);

/*write to file routine to check for input matrix and vector. Uncomment to produce files. Increases considerably overall time and matrix and vector file might be big depending on the precision and number of operation used*/
    // bertini::Mat<bertini::complex> Q=A_b2.fullPivHouseholderQr().matrixQ().eval();//get orthonormal matrix

    // fpA_b1 = fopen ("matrixA_b1.txt" , "a+");//write to file A and b for b1
		// printMat_Matlab_mp(fpA_b1, sprec, A_b1);
    // fclose(fpA_b1);
    // fpb_b1 = fopen ("vectorb_b1.txt" , "a+");
		// printVec_Matlab_mp(fpb_b1, sprec, b_b1);
    // fclose(fpb_b1);

    // std::stringstream buffer;//write to file A and b for b2
    // buffer << A_b2.real().format(MatlabFmtInline)<<"+i*"<<A_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpA_b2 = fopen ("matrixA_b2.txt" , "a+");
    // fprintf(fpA_b2,"%s\n",buffer.str().c_str());
    // buffer.str( std::string() );
    // buffer.clear();

    // fpA_b2.open("matrixA_b2.txt" , std::ios::out | std::ios::app);
    // fpA_b2 <<"\n"<<i<<") "<< A_b2.real().format(MatlabFmtInline)<<"+i*"<<A_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpA_b2.close();
    // fpb_b2.open("vectorb_b2.txt" , std::ios::out | std::ios::app);
    // fpb_b2 <<"\n"<<i<<") "<< b_b2.real().format(MatlabFmtInline)<<"+i*"<<b_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpb_b2.close();

    // buffer <<"\n"<<i<<") "<< b_b2.real().format(MatlabFmtInline)<<"+i*"<<b_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpb_b2 = fopen ("vectorb_b2.txt" , "a+");
    // fprintf(fpb_b2,"%s\n",buffer.str().c_str());//need to change to c++ stream
    // buffer.str("");
    // buffer.clear();
    //
    // buffer << Q.real().format(MatlabFmtInline)<<"+i*"<<Q.imag().format(MatlabFmtInline) << std::endl;
    // fpA_b2 = fopen ("matrixA_b2.txt" , "a+");
    // fprintf(fpA_b2,"%s\n",buffer.str().c_str());
    // buffer.str("");
    // buffer.clear();
/*computation routine for b1. Uncomment the desired solver.*/
    ct1_b1=std::chrono::high_resolution_clock::now();//start timing of one iteration for b1

    // matrixSolve_mp( x_b1, A_b1, b_b1);//b1 solver
    // matrixSolve_cond_num_norms_mp( x_b1,  A_b1,  b_b1,  &cond_num,  &norm_A,  &norm_A_inv);
    matrixSolve_LU_mp( x_b1,  A_b1,  b_b1,  tol,  largeChange);
    // matrixSolve_svd_Least_Squares_mp( x_b1,  A_b1,  b_b1,  tol,  largeChange,  &cond_num,  &norm_A,  &norm_A_inv);  //least square svd
    // matrixSolve_svd_Least_Squares_mp2( x_b1,  A_b1,  b_b1,  mpftol,  mpflargeChange,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_Least_Squares_mp( x_b1,  A_b1,  b_b1,  tol,  largeChange,  &cond_num,  &norm_A,  &norm_A_inv);  // Least square QR
    // matrixSolve_Least_Squares_mp2( x_b1,  A_b1,  b_b1,  mpftol,  mpflargeChange,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_Hessenberg_Least_Squares_mp( x_b1,  A_b1,  b_b1, mpftol, mpflargeChange);
/*computation routine for b2. Uncomment the desired solver.*/
    ct2_b1=std::chrono::high_resolution_clock::now();//end timing of one iteration for b1
    elapsed_seconds_b1 += ct2_b1-ct1_b1;//update total solving timing for b1

    ct1_b2=std::chrono::high_resolution_clock::now();//start timing of one iteration for b2

    // x_b2 = A_b2.colPivHouseholderQr().solve(b_b2);//b2solver
    // x_b2 = A_b2.householderQr().solve(b_b2);
    // x_b2 = A_b2.partialPivLu().solve(b_b2);
    // x_b2 = A_b2.jacobiSvd().solve(b_b2);
    x_b2 = A_b2.partialPivLu().solve(b_b2);

    ct2_b2=std::chrono::high_resolution_clock::now();//end timing of one iteration for b2
    elapsed_seconds_b2 += ct2_b2-ct1_b2;//update total solving timing for b2

    /*write to file routine to check for  result vector. Uncomment to produce files. Increases considerably overall time and matrix and vector file might be big depending on the precision and number of operation used*/
    // fpx_b1 = fopen ("vectorx_b1.txt" , "a+");//write to file solution x for b1
		// printVec_Matlab_mp(fpx_b1, sprec, x_b1);
    // fclose(fpx_b1);
    //
    // fpx_b2.open("vectorx_b2.txt" , std::ios::out | std::ios::app);
    // fpx_b2 <<"\n"<<i<<") "<< x_b2.real().format(MatlabFmtInline)<<"+i*"<<x_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpx_b2.close();

    // buffer << x_b2.real().format(MatlabFmtInline)<<"+i*"<<x_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpx_b2 = fopen ("vectorx_b2.txt" , "a+");
    // fprintf(fpx_b2,"%s\n",buffer.str().c_str());
    // buffer.str("");
    // buffer.clear();

/*routine to output percentage of completion and timing*/
    //Print each percent of loop
      if ( tmp != (int)100*(i+1)/iter )
        {
     	 	tmp = (int)100*(i+1)/iter;
     	 	if ((tmp%1)  == 0)
         {
          std::cout << "Percentage complete: " << tmp<< "%\n";
          std::cout << "elapsed time for linear algebra operations in b1: " << elapsed_seconds_b1.count() << "s\n";//To screen output /*TODO change to write on file providing relevant parameters as well as timing*/
          std::cout << "elapsed time for linear algebra operations in b2: " << elapsed_seconds_b2.count() << "s\n"; }
	       }
    }
    for_end=std::chrono::high_resolution_clock::now();//end total timing of for-loop
    for_laps=for_end-for_start;//calculate total for-loop time

    std::cout << "Percentage complete: 100%\n";
    std::cout << "elapsed time for linear algebra operations in b1: " << elapsed_seconds_b1.count() << "s\n";//To screen output /*TODO change to write on file providing relevant parameters as well as timing*/
    std::cout << "elapsed time for linear algebra operations in b2: " << elapsed_seconds_b2.count() << "s\n";
    std::cout << "elapsed time during for-loop: " << for_laps.count() << "s\n";

    clear_mat_mp(A_b1);//free memory for b1 vectors and matrix /*Not needed in Eigen?*/
    clear_vec_mp(x_b1);
    clear_vec_mp(b_b1);

    mpf_clear(mpftol);//clear mpfr parameters
    mpf_clear(mpflargeChange);

    clearMP();// clears bertini 1 globals
    return 0;

}
