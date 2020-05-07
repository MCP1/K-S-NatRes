EQUATION("_teste2")
v[0]=3;
RESULT( v[0] )

EQUATION( "_Ataup" )
/*
Productivity of the new vintage of machines when employed for production.
Also updates '_Btau'.
*/
double Btau = VL( "_Btaup", 1 );					// previous period productivity
double xi = VS( PARENT, "xi" );
// normalized workers on R&D of the firm
double L1rdN = VL( "_L1rdp", 1 ) * VS( LABSUPL2, "Ls0" ) / VLS( LABSUPL2, "Ls", 1 );
// innovation process (success probability)
v[1] = 1 - exp( - VS( PARENT, "zeta1")*xi);//*L1rdN); THERE ARE PROBLEMS MULTIPLYING BY L1rndN
// normalized workers on R&D of the firm
if ( bernoulli( v[1] ) )						// innovation succeeded?
{
	double x1inf = VS( PARENT, "x1inf" );		// lower beta inno. draw support 
	double x1sup = VS( PARENT, "x1sup" );		// upper beta inno. draw support 
	double alpha1 = VS( PARENT, "alpha1" );		// beta distrib. alpha parameter
	double beta1 = VS( PARENT, "beta1" );		// beta distrib. beta parameter
		// new final productivity (A) from innovation
	v[2] = CURRENT * ( 1 + x1inf + beta( alpha1, beta1 ) * ( x1sup - x1inf ) );
		// new production productivity (B) from innovation
	v[3] = Btau * ( 1 + x1inf + beta( alpha1, beta1 ) * ( x1sup - x1inf ) );
}
else
	v[2] = v[3] = 0;							// innovation failure


//Add Immitation Module
// Need to completely adapt it to a public firm - There are problems here

v[4] = 1 - exp( - VS( PARENT, "zeta2" ) * ( 1 - xi ) );//*L1rdN 

if ( bernoulli( v[4] ) )						// imitation succeeded?
{
	k = VS( PARENT, "F1" );					//CHECK THIS	// number of firms in sector 1 - See if I don't have to add the public firms to F1 as well
	dblVecT imiProb( k );						// vector for tech distance - CHECK THIS FUNCTION.
	
	v[5] = i = 0;								// inverse distance/firm accum.
	CYCLES( PARENT, cur, "Firm1p" )				// 1st run: abs. inv. distance
		if ( cur == p )
			imiProb[ i++ ] = 0;					// can't self-imitate
		else
		{
			v[6] = sqrt( pow( VLS( cur, "_Btaup", 1 ) - Btau, 2 ) +
						 pow( VLS( cur, "_Ataup", 1 ) - CURRENT, 2 ) );
			v[5] += imiProb[ i++ ] = ( v[6] > 0 ) ? 1 / v[6] : 0;
		}

	if ( v[5] > 0 )
	{
		v[7] = i = 0;							// probabilities/firm accum.
		CYCLES( PARENT, cur, "Firm1p" )			// 2nd run: cumulative imi. prob.
		{
			v[7] += imiProb[ i ] / v[5];		// normalize to add up to 1
			imiProb[ i++ ] = v[7];
		}
			
		// draw a firm to imitate according to the distance probabilities
		j = upper_bound( imiProb.begin( ), imiProb.end( ), RND ) - imiProb.begin( );
		
		if ( j < k )							// at least one firm reachable?
		{
			cur = TSEARCHS( PARENT, "Firm1p", j + 1 );// get pointer to firm
			
			v[8] = VLS( cur, "_Ataup", 1 );		// get imitated firm productivities
			v[9] = VLS( cur, "_Btaup", 1 );
		}
		else
			v[8] = v[9] = 0;					// imitation failure
	}
	else
		v[8] = v[9] = 0;						// imitation failure
}
else
	v[8] = v[9] = 0;		

// End of Immitation Module


// select best option between the three options (current/innovation/imitation)
v[0] = CURRENT;									// current technology
v[10] = Btau;
v[11] = v[12] = 0;
if (v[2] * v[3] > v[0] * v[10] )				// is innovation better?
{
	v[0] = v[2];								// new Atau
	v[10] = v[3];								// new Btau
	v[11] = 1;									// innovation succeeded
	v[12] = 0;									// no imitation
}
if ( v[8] * v[9] > v[0] * v[10] )				// is imitation better (yet)?
{
	v[0] = v[8];
	v[10] = v[9];
	v[11] = 0;									// no innovation
	v[12] = 1;									// imitation succeeded
}
WRITE( "_Btaup", v[10] );
WRITE( "_innp", v[11] );
WRITE( "_imip", v[12] );
RESULT( v[0] )
// I Need _Btaup _L1rdp _inip _imip

EQUATION_DUMMY( "_Btaup", "" )
/*
Productivity of labor in producing the new vintage of machines
Updated in '_Atau'
*/

EQUATION_DUMMY( "_imip", "" )
/*
Imitation success (1) or failure (0) for firm in capital-good sector
Updated in '_Atau'
*/

EQUATION_DUMMY( "_innp", "" )
/*
Innovation success (1) or failure (0) for firm in capital-good sector
Updated in '_Atau'
*/


EQUATION( "_L1rdp" )
/*
R&D labor employed by firm in capital-good sector
*/
v[0]= V("_L1dRDp") * VS( PARENT, "L1rd" ) / VS( PARENT, "L1dRD" );
RESULT( v[0] )
// I Need _L1dRDp

EQUATION( "_L1dRDp" )
/*
R&D labor demand of firm in capital-good sector
*/
RESULT( V("_RDp") / VS( LABSUPL2, "w" ) )
//I Need _RD

EQUATION( "_RDp" )
/*
R&D expenditure of firm in capital-good sector
*/
v[1] = VL( "_S1p", 1 );					// sales in previous period
v[2] = VS( PARENT, "nu" );						// R&D share of sales
if ( v[1] > 0 )
// This is where I add Francesco's first experiment
	v[0] = v[2] * v[1];
else											// no sales
	// keep current expenditure or a share of available cash
	v[0] = CURRENT;
RESULT( v[0] )
//I Need _S1p

EQUATION( "_S1p" )
/*
Sales of firm in capital-good sector
*/
RESULT( V( "_p1p" ) * V( "_Q1ep" ) )
// I Need _p1 and _Q1e

// To _p1p
EQUATION( "_p1p" )
/* Price of good of firm in capital-good sector */
RESULT( ( 1 + VS( PARENT, "mu1" ) ) * V( "_c1p" ) )
// I Need _c1p

EQUATION( "_c1p" )
/* Unit cost of production (1 machine) of firm in capital-good sector.
Use sectoral average (pool of sector 1 workers). */
V("_Ataup");									// ensure innovation process ok
RESULT( VS( LABSUPL2, "w" ) / ( V( "_Btaup" ) * VS( PARENT, "m1" ) ) )

// To _Q1e
EQUATION( "_Q1ep" )
/* Effective output of firm in capital-good sector */
v[0] = V( "_Q1p" );								// planned production
v[1] = V( "_L1p" );								// effective labor available
v[2] = V( "_L1dp" );								// desired total workers
if ( v[1] >= v[2] )
	END_EQUATION( v[0] );						// produce as planned
v[3] = V( "_L1rdp" );							// effective R&D workers
v[4] = V( "_L1dRDp" );							// desired R&D workers
v[5] = 1 - ( v[1] - v[3] ) / ( v[2] - v[4] );	// adjustment factor
// adjust all pending orders, supplying at least one machine
CYCLE( cur, "Clip" )
	if ( VS( cur, "_tOrdp" ) == T )				// order in this period?
	{
		v[6] = VS( cur, "_nOrdp" ) - VS( cur, "_nCanp" );// existing net orders
		v[0] -= v[7] = min( floor( v[6] * v[5] ), v[6] - 1 );
		INCRS( cur, "_nCanp", v[7] );
	}
RESULT( max( v[0], 0 ) )						// avoid negative in part. cases
// I Need _Q1p _L1p _L1dp

// First _L1dp
EQUATION( "_L1dp" )
/*
Labor demand of firm in capital-good sector
Includes R&D labor
*/
RESULT( V( "_L1dRDp" ) + V( "_Q1p" )  / ( V( "_Btaup" ) * VS( PARENT, "m1" ) ) )

// Second _L1p
EQUATION( "_L1p" )
/*
Labor employed by firm in capital-good sector Includes R&D labor
*/
v[1] = VS( PARENT, "L1" );						// total labor in sector 1
v[2] = VS( PARENT, "L1d" );						// total labor demand in sector 1
v[3] = VS( PARENT, "L1rd" );					// R&D labor in sector 1
v[4] = VS( PARENT, "L1dRD" );					// R&D labor demand in sector 1
v[5] = V( "_L1dp" );								// firm total labor demand
v[6] = V( "_L1rdp" );							// firm R&D labor
v[7] = V( "_L1dRDp" );							// firm R&D labor demand
// total workers in firm after possible shortage of workers
RESULT( min( v[5], v[6] + ( v[2] > v[4] ? ( v[5] - v[7] ) * ( v[1] - v[3] ) / 
										  ( v[2] - v[4] ) : 0 ) ) )
										  
// Third, _Q1p
EQUATION( "_Q1p" )
/*
Planed production for a firm in capital-good sector
*/
v[1] = V( "_D1p" );								// potential production (orders)
v[2] = V( "_cred1p" );							// available credit
v[3] = VL( "_NW1p", 1 );							// net worth (cash available)
v[4] = V( "_c1p" );								// unit cost
v[5] = V( "_p1p" );								// machine price
v[6] = V( "_RDp" );								// R&D costs still to pay
v[7] = v[1] * ( v[4] - v[5] ) + v[6];			// cash flow to fulfill orders

if ( v[7] < 0 || v[7] <= v[3] - 1 )				// firm can self-finance?
{
	v[0] = v[1];								// plan the desired output
	v[3] -= v[7];								// remove flow from cash
}
else
{
	if ( v[7] <= v[3] - 1 + v[2] )				// possible to finance all?
	{
		v[0] = v[1];							// plan the desired output
		v[8] = v[9] = v[7] - v[3] + 1;			// finance the difference
		v[3] = 1;								// keep minimum cash
	}
	else										// credit constrained firm
	{
		// produce as much as the available finance allows
		v[0] = floor( ( v[3] - 1 - v[7] + v[2] ) / v[4] );// max possible
		v[0] = min( max( v[0], 0 ), v[1] );		// positive but up to D1
		
		v[8] = v[2];							// finance what is possible
		v[9] = v[6] - v[3] + 1;					// desired credit
	
		if ( v[0] == 0 )
		{
			v[10] = 1;							// all orders canceled
			v[3] -= v[6] - v[2];				// let negative NW (bankruptcy)
		}
		else
		{
			v[10] = 1 - v[0] / v[1];			// machine shortage factor
			v[3] = 1;							// keep minimum cash
		}
				// shrink or cancel all exceeding orders
		CYCLE( cur, "Clip" )
			if ( VS( cur, "_tOrdp" ) == T )		// order in this period?
			{
				if ( v[10] == 1 )				// bankruptcy?
					INCRS( cur, "_nCanp", VS( cur, "_nOrdp" ) );
				else
					INCRS( cur, "_nCanp", floor( VS( cur, "_nOrdp" ) * v[10] ) );
			}
	}
		// TO DO TO DO TO DO TO DO TO DO on support.h
		// TO DO TO DO TO DO TO DO TO DO
		//update_debt1p( p, v[9], v[8] );				// update firm debt
		// TO DO TO DO TO DO TO DO TO DO
		// TO DO TO DO TO DO TO DO TO DO
}
// provision for revenues and expenses
WRITE( "_NW1p", v[3] );							// update the firm net worth
RESULT( v[0] )

// I have to add _D1p _cred1p _NW1p 
  
// First NW1p
EQUATION_DUMMY( "_NW1p", "" )
/*
Net wealth (free cash) of firm in capital-good sector
Updated in '_Q1' and '_Tax1'
*/

// Second _D1p
EQUATION( "_D1p" )
/*
Potential demand (orders) received by a firm in capital-good sector
*/
VS( CONSECL2, "Id" );							// make sure all orders are sent CHECK WHAT ARE THESE IDs of consumption good firms
j = v[0] = 0;									// machine/active customer count
CYCLE( cur, "Clip" )
{
	if ( VS( cur, "_tOrdp" ) == T )				// order in this period?
	{
		v[0] += VS( cur, "_nOrdp" );
		++j;
	}
}
WRITE( "_BCp", j );
RESULT( v[0] )

// Dummy Equation _BCp
EQUATION_DUMMY( "_BCp", "" )
/*
Number of buying clients for firm in capital-good sector
Updated in '_D1'
*/

// Third, _cred1p
EQUATION( "_cred1p" )
/*
Bank credit available (new debt) to firm in capital-good sector
Function called multiple times in single time step
*/
v[1] = V( "_Deb1p" );							// current firm debt
v[2] = V( "_Deb1maxp" );							// maximum prudential credit

/*  CHECK THIS PART, THE BANKS ASSOCIATED TO THE FIRMS. A PROBLEM FOR ME UNTIL NOW
if ( v[2] > v[1] )								// more credit possible?
{
	v[0] = v[2] - v[1];							// potential free credit
	
	cur = HOOK( BANK );							// firm's bank - What does this mean?
	v[3] = VS( cur, "_TC1free" );				// bank's available credit // CHECK THIS PART!!!!!!
	
	if ( v[3] > -0.1 )							// credit limit active
		v[0] = min( v[0], v[3] );				// take just what is possible
}
else 
*/
	v[0] = 0;									// no credit available

RESULT( v[0] )
// colocar _Deb1 e _Deb1max

// First _Deb1max
EQUATION( "_Deb1maxp" )
/*
Prudential maximum bank debt of firm in capital-good sector
*/
// maximum debt allowed to firm, considering net worth and operating margin
v[5] = VS( FINSECL2, "Lambda" ) * max( VL( "_NW1p", 1 ), 
									   VL( "_S1p", 1 ) - VL( "_W1p", 1 ) );   
// apply an absolute floor to maximum debt prudential limit
v[0] = max( v[5], VS( FINSECL2, "Lambda0" ) * VLS( PARENT, "PPI", 1 ) / 
				  VS( PARENT, "PPI0" ) );
WRITE( "_cred1cp", 0 );							// reset constraint for period
RESULT( v[0] )
// Do _W1p and dummy equation _cred1cp

EQUATION_DUMMY( "_cred1cp", "" )
/*
Credit constraint for firm in capital-good sector
Updated in '_Deb1max', '_Q1'
*/

// Now W1p
EQUATION( "_W1p" )
/*
Total wages paid by firm in capital-good sector
*/
RESULT( V( "_L1p" ) * VS( LABSUPL2, "w" ) )

//Part of Deb1p
EQUATION_DUMMY( "_Deb1p", "" )
/*
Stock of bank debt of firm in capital-good sector
Updated in '_Q1', '_Tax1'
*/

EQUATION( "_NCp" )
/*
Number of new client firms in the period.
Also creates the client-supplier connecting objects.
*/

firmMapT firms = V_EXTS( GRANDPARENT, country, firm2map );// list with all firms
h = firms.size( );								// number of firms in sector 2
k = V( "_HCp" );									// number of historical clients
j = V( "_ID1p" );								// current firm ID

CYCLE( cur, "Clip" )								// remove historical clients
	firms.erase( ( int ) VS( cur, "_IDcp" ) );
	
i = ceil( VS( PARENT, "gamma" ) * k );			// new clients in period
i = v[0] = min( max( i, 1 ), min( h - k, firms.size( ) ) );// between [1, F2 - HC]

v[1] = max( 1, ceil( h / VS( PARENT, "F1" ) ) );// firm fair share

if ( k + i < v[1] )								// ensure at least fair share
	i = v[0] = v[1] - k;

// build vector of all target firms (not yet clients)

//**************************SOMETHING NEEDS TO BE DONE IN THIS PART OTHERWISE THE MODEL DOES NOT RUN***********************/

//vector < firmPairT > targets( firms.begin( ), firms.end( ) );
/*
// draw new clients from target list, updating the list after each draw
for ( ; i > 0; --i )
{
	h = uniform_int( 0, targets.size( ) - 1 );	// draw a index in the list
	firmPairT client = targets[ h ];			// retrieve drawn map pair
	targets.erase( targets.begin( ) + h );		// remove drawn firm from list

	// create the brochure/client interconnected objects
	send_brochure( j, p, client.first, client.second );
}
*/
RESULT( v[0] )

// Need to create _HC
EQUATION( "_HCp" )
/*
Number of historical client firms from consumer-good sector.
Also removes old, non-buying clients.
*/
i = 0;											// client counter
CYCLE_SAFE( cur, "Clip" )						// remove old clients
{
	if ( VS( cur, "_tSelp" ) < T - 1 )			// last selection is old?
	{
	//THERE ARE PROBLEMS WITH THIS POINTER WHEN I RUN THE MODEL
	//	DELETE( SHOOKS( cur ) );				// remove supplier brochure entry
	//	DELETE( cur );							// remove client entry
	}
	else
		++i;
}
RESULT( i )


// Still missing Tax1 J01 pi1 f1 div1 qc1,cscores


EQUATION( "_Tax1p" )
/*
Total tax paid by firm in capital-good sector
Also updates final net wealth on period
*/
v[1] = V( "_Pi1p" );								// firm profit in period
v[2] = VS( GRANDPARENT, "tr" );					// tax rate
if ( v[1] > 0 )									// profits?
{
	v[0] = v[1] * v[2];							// tax to government
	v[4] = VS( PARENT, "d1" ) * ( v[1] - v[0] );// dividend to shareholders
}
else
	v[0] = v[4] = 0;							// no tax/dividend on losses
WRITE( "_Div1p", v[4] );							// save period dividends

// compute free cash flow
v[6] = v[1] - v[0] - v[4];

// remove from net wealth the provision for revenues and expenses
v[7] = INCR( "_NW1p", V( "_Q1p" ) * ( V( "_c1p" ) - V( "_p1p" ) ) + V( "_RDp" ) );
if ( v[6] < 0 )									// must finance losses?
{
	if ( v[7] >= - v[6] + 1 )					// can cover losses with reserves?
		INCR( "_NW1p", v[6] );					// draw from net wealth
	else
	{
		v[8] = V( "_cred1p" );					// available credit
		v[9] = - v[6] - v[7] + 1;				// desired finance
		
		if ( v[8] >= v[9] )						// can finance losses?
		{
			update_debt1( p, v[9], v[9] );		// finance all
			WRITE( "_NW1p", 1 );					// minimum net wealth
		}
		else
		{
			update_debt1( p, v[8], v[8] );		// take what is possible // TO CHECK
			INCR( "_NW1p", v[6] - v[8] );		// let negative NW (bankruptcy exit)
		}					
	}
}
else											// pay debt with available cash
{
	v[10] = V( "_Deb1p" );						// current debt
	
	if ( v[10] > 0 )							// has debt?
	{
		if ( v[6] > v[10] )						// can repay all debt and more
		{
			update_debt1( p, 0, - v[10] );		// zero debt
			INCR( "_NW1p", v[6] - v[10] );		// save the rest
		}
		else
			update_debt1( p, 0, - v[6] );		// repay part of debt //CHECK
	}
	else
		INCR( "_NW1p", v[6] );					// save all
}	
RESULT( v[0] )

// Needs _Div1p _pi1p

EQUATION_DUMMY( "_Div1p", "" )
/*
Dividends paid by firm in capital-good sector
Updated in '_Tax1'
*/

EQUATION( "_Pi1p" )
/* Profit of firm in capital-good sector */
v[1] = V( "_S1p" ) - V( "_W1p" );					// gross operating margin
v[2] = VS( FINSECL2, "rD" ) * VL( "_NW1p", 1 );	// financial income
// firm effective interest rate on debt
v[3] = VLS( FINSECL2, "rDeb", 1 ) * ( 1 + ( VL( "_qc1p", 1 ) - 1 ) * VS( FINSECL2, "kConst" ) ); 
v[4] = v[3] * VL( "_Deb1p", 1 );					// interest to pay
RESULT( v[1] + v[2] - v[4] )					// firm profits before taxes

//Falta _qc1p

EQUATION_DUMMY( "_qc1p", "cScores" )
/*
Credit class of firm in sector 1 (1,2,3,4)
Updated in 'cScores'
*/

EQUATION( "_JO1p" )
/*
Open job positions for a firm in capital-good sector
*/
RESULT( max( V( "_L1dp" ) - VL( "_L1p", 1 ), 0 ) )

EQUATION( "_f1p" )
/*
Market share of firm in capital-good sector
Keep market share if sector didn't produce
*/
v[1] = VS( PARENT, "Q1e" );
RESULT( v[1] > 0 ? V( "_Q1ep" ) / v[1] : CURRENT )







// To Do
// Update_debt1 - Adapt it to the public firm on support.h. I have to understand how it works. 
