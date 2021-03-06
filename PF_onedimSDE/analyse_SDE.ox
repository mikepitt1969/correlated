#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>

Parzen( const lag);
var_inflate(const lag, const mX);

/* Analyses output from the main file which is called:
  pmcmc_SQRT_nolev.ox OR  pmcmc_SQRT_lev.ox
*/

main()
{

   
   // columns are: difference in logL (prop - current), logL, accept rate, 4 parameters (mu, phi = exp(-nu), omega, squiggly X (leverage)
   // 7 cols
	
	
	decl mX =  loadmat( "output_nolev.mat" )[][]; 			 
	decl r = rows(mX), c =  columns(mX);
	print(r~c);	  //, mX[][1:10]
	decl RUN_g = 300, BURN_IN = 3000;
		
    Draw(0,mX[3:6][0:(c-1)] , 0, 1);	//	 i+3  	for (decl i =0; i < 4; i++) 

	DrawTitle(0, " CPM: Draws of parameters" );	// DrawTitle(3, " R=Z'-Z" );


	decl lag = 820;	   //c=10000;
	decl acfX = acf(mX[3:5][3500:(c-1)]', lag)[0:lag][]';
	Draw(1, acfX[0:2][], 0, 1);	 DrawTitle(1, " CPM: Correlogram for mu, phi, omega " );
	 
	print(1 + 2 * sumc(acfX'));
	//lag = 120;
	acfX = acf(mX[6][3500:(c-1)]', lag)[0:lag][]';
	Draw(2, acfX[0][], 0, 1);	 DrawTitle(2, " CPM: Correlogram for leverage parameter " );
	
	Draw(3, mX[0][(RUN_g+1):(BURN_IN-1)] , 0, 1);   DrawTitle(3, " logL proposed minus current value " );
	Draw(4,  mX[1][10:(c-1)], 0, 1);   DrawTitle(4, " logL " );
    Draw(5, mX[2][10:(c-1)] , 0, 1);   DrawTitle(5, " Av acceptance prob " );
	acfX = acf(mX[0][(RUN_g+1):(BURN_IN-1)]', lag)[0:lag][]';
	Draw(6, acfX[][], 0, 1);	 DrawTitle(6, " Correlogram for diff in logLik (at fixed par) " );
	DrawDensity(7, mX[0][(RUN_g+1):(BURN_IN-1)]," W =Z'-Z",  1, 1, 1);
	   	ShowDrawWindow();


		
	decl sigma2= varr(	mX[1][0:RUN_g-1] );
	print( sigma2 );
		exit(1);


		
	print(1 + 2 * sumc(acfX'));

//   Draw(2, mX[2][0:99], 0, 1);
	DrawDensity(7, mX[1][10:(c-1)]," Z under g()",  1, 1, 1);
	DrawDensity(2, mX[2][100:(c-1)]," Z under pi()",  3, 3, 3);

	Draw(3, mX[5][4000:(c-1)], 0, 1);   DrawTitle(3, " R=Z'-Z" );
	
	DrawDensity(4, mX[5][4000:(c-1)]," Histogram R=Z'-Z ",  1, 1, 1);
	lag = 100;
		acfX = acf(mX[5][2000:(c-1)]', lag)[0:lag][]';
	Draw(5, acfX[0][], 0, 1);	 DrawTitle(5, " correlogram of R=Z'-Z" );

	
	ShowDrawWindow();
	SaveDrawWindow("a.eps")	;
}




	//	Draw(1,mX0[0][] , 0, 1);
//	DrawDensity(2, mX[0][100:(c-1)],"",  1, 1, 1);
//	DrawDensity(3, mX0[0][],"",  1, 1, 1);
   	//DrawDensity(0, mX[0][100:(c-1)],"",  1, 1, 1);
//	for (decl j=0; j < 7; j++)	Draw(j+6, mX[j][5:2000],0,1); // Draw(1, mX[0:6][], 0, 1)	;
//	print(meanr(mX[6][5000:(c-1)])~varr(mX[6][5000:(c-1)]) );
//	decl kappa = sqrt( varr(mX[5][5000:(c-1)]) );
//	print( 2 * probn(-kappa/2) );
//	DrawDensity(7,mX[5][5000:(c-1)] ,"",  1, 1, 1);
//
//	 decl acfX = acf(mX[0][5100:(c-1)]', lag)[0:lag][]';
//	 Draw(8, acfX[0][], 0, 1);
//	decl mC =  loadmat( "sojourns_med_d2.mat" )[][]';
//	 r = rows(mC), c =  columns(mC);
//	print(mC[0][0:10]);
//	Draw(6, mC[1][2000:(c-1)], 0, 1);  DrawDensity(7,mC[1][2000:(c-1)],"",  1, 1, 1);
   	//	 print(	"accept ",  (mX[2][c-1] * c -40)/(c-40)~
	//	 2 * probn(-kappa/2)~ 2 * probn( -sqrt( varr(mX[0][0:(c-1)]))/2) );
	//decl acfY = acf(mX[0][100:(c-1)]', 100)[1:100][]';
	// Draw(2, acfY[0][], 0, 1);
	//Draw(1, mX[1][10:(c-1)], 0, 1)	;
	//Draw(2, mX[2][0:(c-1)], 0, 1)	;
	//print(mX[][0:10]');
//	decl lag = 700;
//	decl acfX = acf(mX[3:6][5100:(c-1)]', lag)[0:lag][]';
//	Draw(2, mX[3][], 0, 1)	; Draw(3, acfX[0][], 0, 1);//DrawCorrelogram(3,mX[3][1000:(c-1)] , "", 2000);	
//	Draw(4, mX[4][], 0, 1)	; Draw(5, acfX[1][], 0, 1);//DrawCorrelogram(5,mX[4][1000:(c-1)] , "", 2000);	
//	Draw(6, mX[5][], 0, 1)	;  Draw(7, acfX[2][], 0, 1);//DrawCorrelogram(7,mX[5][1000:(c-1)] , "", 2000);	
//	Draw(8, mX[6][], 0, 1)	;  Draw(9, acfX[3][], 0, 1);//DrawCorrelogram(9,mX[6][3000:(c-1)] , "", 2000);	
//	ShowDrawWindow();
//
//	decl index = range(1, c-201);
//	print( var_inflate(700, mX[3:6][1100:(c-1)]) ) ;
//	//print(  var_inflate(1000, mX[6][3000:(c-1)]) );

	
	/*
	//decl vx = range(1, 1500) ./100 -10;
	//print(r~c);
	decl kappa1 = sqrt( varr(mX[0][1:299]) ) ;
	print(kappa1~-2 *meanr(mX[0][1:299]));	  decl kappa2 =  sqrt(varr(mX[1][201:(c-1)]) );
	//Draw(4, (mX[2][201:(c-1)] .*(200+index) -200)./index | 2*probn(-kappa2/sqrt(2)) .* ones(1,c-201), 0, 1)	;
	ShowDrawWindow();
	print(meanr( mX[3:6][6000:(c-1)])~sqrt( varr( mX[3:6][6000:(c-1)]) ) );
	//	  //-2 *meanr(mX[0][2000:(c-1)])
	//print(kappa1~kappa2~ 2*probn(-kappa1/sqrt(2))~2*probn(-kappa2/sqrt(2)));
	*/	 
	

 	//Draw(0, mX0[1][0:1500], 0, 1)	;
	//DrawDensity(0,mX0[1][] ,"",  1, 1, 1);
	//Draw(1, mX[1][0:(c-1)], 0, 1)	;
	//DrawDensity(1,mX[1][1500:(c-1)] ,"",  1, 1, 1);

	//Draw(2, mX[0][0:(c-1)], 0, 1)	;
	//DrawDensity(3, mX[0][0:(c-1)],"",  1, 1, 1);
	// decl acfY = acf(mX[0][1500:(c-1)]', 50)[0:50][]';
	//Draw(4, acfY[0][], 0, 1);
	//decl index = range(1, c-41);
	//decl kappa =  sqrt(-2 * meanr(mX[0][0:(c-1)]) );
//Draw(5, ( mX[2][41:(c-1)].*(41+index) -40 )./index | 2 * probn(-kappa/2) .* ones(1, c-41)
	//     | 2 * probn( -sqrt( varr(mX[0][0:(c-1)]))/2) .* ones(1, c-41), 0, 1)	;
	//Draw(4, 2 * probn(-kappa/2) .* ones(1, c-41), 0, 1)	;

Parzen( const lag)
{
  decl working, i;
  decl parzen = zeros(1, lag);
  for (i = 0; i < lag; i++)
  {
      working = double(i)/double(lag);
      if (working < 0.5)		  
		parzen[][i] = 1.0 - (6.0 * (working)^2) + (6.0 * (working)^3);
      else
		parzen[][i] = 2.0 * (1.0 - working)^3 ;
  }
  return parzen;
}

// mX 1*T.
var_inflate(const lag, const mX)
{
  decl acfX = acf(mX', lag)[1:lag][];
  decl par = ones(1, lag);//Parzen(lag);
  return 1 + 2 * sumc(acfX .* par');
}