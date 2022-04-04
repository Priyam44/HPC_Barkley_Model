/**
 * Author   :  Priyam Gupta
 * Username :  pg221
 * CID#     :  02110124
 * For      :  AERO70011 High Perfomance Computing Coursework
 * Objective : Parallelised implementation using OpenMP for solving
               Barkley Reaction Diffusion model uisng Forward Euler
               in time and Central difference in space schemes
 **/

class ReactionDiffusion{

  //Data Members
  int Nx;            //Number of Grid points in x
  int Ny;            //Number of Grid points in y
  int Nx_intvl;      //Number of Grid intervals in x
  int Ny_intvl;      //Number of Grid intervals in y
  int bNx;           //Extended number of Grid points to utilise normal scheme for the Boundary Condition in x
  int bNy;           //Extended number of Grid points to utilise normal scheme for the Boundary Condition in y
  int nthreads;      //Number of Threads to be spawned

  double dt;         //Timestep
  double T;          //Total integration time
  double a;          //Parameter a
  double b;          //Parameter b
  double eps;        //Paramter epsilon
  double mu1;        //Parameter mu1
  double mu2;        //Paramter mu2
  double dx;         //Spatial Element Size in x
  double dy;         //Spatial Element Size in y
  double Lx;         //Length of Domain in x
  double Ly;         //Length of Domain in y
  double sigma1;     //Discritization Coefficient for u
  double sigma2;     //Discritization Coefficienti for v
  double  f_1  ;     //Reaction Element for u
  double  f_2  ;     //Reaction Element for v

  double * un  ;     //u matrix
  double * vn  ;     //v matrix
  double * un1 ;     //u matrix at subsequent timestep
  double * vn1 ;     //v matrix at subsequent timestep


  public:
  //////////////////////////////////////////////////////////////////////////////
  /*Assigns the required Parameter values to corresponding
    ReactionDiffusion Class data members
   *Takes all the required Parameters as input
   *Calculates the variables dependent on input parameters
    required for solving the partial differential equation
  */
  void SetParameters(int Nx, int Ny, double a, double b, double eps, double mu1, double mu2, double T, double dt, double dx, double dy, int nthreads){
    this-> Nx  = Nx;
    this-> Ny  = Ny;
    this-> a   = a;
    this-> b   = b;
    this-> eps = eps;
    this-> mu1 = mu1;
    this-> mu2 = mu2;
    this-> T   = T;
    this-> dt  = dt;
    this-> dx  = dx;
    this-> dy  = dy;
    this-> nthreads = nthreads;

    bNx        = Nx + 2;
    bNy        = Ny + 2;
    Nx_intvl   = Nx - 1;
    Ny_intvl   = Ny - 1;
    Lx         = dx*(Nx_intvl);
    Ly         = dy*(Ny_intvl);
    sigma1     = mu1*dt/(dx*dx);
    sigma2     = mu2*dt/(dy*dy);

    un         = new double [bNx * bNy];
    vn         = new double [bNx * bNy];
    un1        = new double [bNx * bNy];
    vn1        = new double [bNx * bNy];

  }
  //////////////////////////////////////////////////////////////////////////////
  /*Intializes u and v matrices to the given Initial Conditions
  */
  void SetInitialConditions(){

    int i,j,ij;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for private(i,ij) collapse(2)
    for (j=0; j<bNx; ++j){
      for (i=0; i<bNy; ++i){
        //Precalculating index
        ij = j*bNy + i;

        //Assigning Initial Conditions to the u and v matrices
        un[ij] =  (i*dy) > Ly*0.5 + 1 ? 1 : 0;
        vn[ij] =  (j*dx) < Lx*0.5 + 1 ? a*0.5 : 0;

      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  void TimeIntegrate(){
    //Time Integrate
    int num_timesteps = T/dt;
    for(int t=0;t<num_timesteps;++t){
      if(t%40000 == 0){
        cout <<"Timestep = " << t <<endl;
      }

      int i,j,ij;                    //Initialzing Loop Vaiables

      //Parallel Loop running over the grid
      omp_set_num_threads(nthreads);
      #pragma omp parallel for private(i,ij,f_1,f_2) collapse(2)
      for (j=1; j < bNx-1; ++j){
        for (i=1; i < bNy-1; ++i){

          //Precalculating Index
          ij = j * bNy + i;
          //Calculating Reaction Term
          f_1 =  eps * un[ij] * (1 - un[ij]) * (un[ij] - (vn[ij] + b)/a);
          f_2 =  un[ij]*un[ij]*un[ij] - vn[ij];

          //Implementing Forward Euler in Time and Central Difference in Space Scheme
          un1[ij] = un[ij] * (1 - 4 * sigma1) + sigma1 * (un[ij + 1] + un[ij - 1] + un[ij + bNy] + un[ij - bNy]) + dt * f_1;
          vn1[ij] = vn[ij] * (1 - 4 * sigma2) + sigma2 * (vn[ij + 1] + vn[ij - 1] + vn[ij + bNy] + vn[ij - bNy]) + dt * f_2;

          //Updating Left Dummy Variables
          if(j == 1){
            un1[ij - bNy] = un1[ij];
            vn1[ij - bNy] = vn1[ij];
          }
          //Updating Right Dummy Variables
          if(j == Nx){
            un1[ij + bNy] = un1[ij];
            vn1[ij + bNy] = vn1[ij];
          }
          //Updating Top Dummy Variables
          if(i == Ny){
            un1[ij + 1] = un1[ij];
            vn1[ij + 1] = vn1[ij];
          }
          //Updating Bottom Dummy Variables
          if(i == 1){
            un1[ij - 1] = un1[ij];
            vn1[ij - 1] = vn1[ij];
          }
        }
      }

      //Parqallel copying of Un1 into Un for next time step calculation
      omp_set_num_threads(nthreads);
      #pragma omp parallel for private(i) collapse(2)
      for(j = 0; j < bNx; ++j){
        for(i = 0; i < bNy; ++i){
          un[j * bNy + i] = un1[j * bNy + i];
          vn[j * bNy + i] = vn1[j * bNy + i];
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  /*Prints Matrix
    Input: Matrix in Column Major Format, Number of Rows, Number of Columns
  */
  void printmatrix(double *A, int r, int c){

          for(int i = 0 ;i<r;++i){
            for (int j=0; j<c;++j){
                  cout<<setw(15)<<A[j*r + i];
          }cout<<endl;
      }
  }

  //////////////////////////////////////////////////////////////////////////////
  //Returns U matrix
  double * Get_u (){
    return un;
  }
  //Returns V matrix
  double * Get_v(){
    return vn;
  }
  //Prints out the important Parameters being used by the program
  void print_para(){
    cout << endl;
    cout << "###################################################" << endl;
    cout << "Nx: " << Nx << endl;
    cout << "Ny: " << Ny << endl;
    cout << "T: " << T << endl;
    cout << "a: " << a << endl;
    cout << "b: " << b << endl;
    cout << "eps: " << eps << endl;
    cout << "mu1: " << mu1 << endl;
    cout << "mu2: " << mu2 << endl;
    cout << "dx: " << dx << endl;
    cout << "dy: " << dy << endl;
    cout << "Nx_intvl: " << Ny_intvl << endl;
    cout << "Lx: " << Lx << endl;
    cout << "Ly: " << Ly << endl;
    cout << "sigma1: " << sigma1 << endl;
    cout << "sigma2: " << sigma2 << endl;
    cout << "###################################################" << endl;
    cout << endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  /*Writes Velocity Data into file
   *Takes desired filename as input
  */
  void write(string filename){

    ofstream fout;
    fout.open(filename);

    for(int i=1; i<bNy-1; ++i){
      for (int j=1; j<bNx-1; ++j){
        fout << j << " " << i << " " << un[j*bNy + i] << " " << vn[j*bNy + i]<< endl;
      }fout << endl;
    }
    fout.close();
  }
  //////////////////////////////////////////////////////////////////////////////

  /*Destructor
   *Deletes the allocated array memory during runtime
  */
};
