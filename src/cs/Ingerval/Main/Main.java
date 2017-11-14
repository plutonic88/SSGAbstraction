package cs.Ingerval.Main;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;
import groupingtargets.ClusterTargets;
import groupingtargets.GroupingTargets;

public class Main {


	/**
	 * 
	 * 1. Assume that S is a 2-D array where a Sudoku grid is stored.
2. sum = S.length * (S.length + 1) / 2;
3. for r =0 to S.length - 1
4. missingNum = sum;
5. for i = 0 to S.length – 1
6. if (S[r][i] contains a known value)
7. missingNum - = S[r][i];
8. print("The missing number in row: " + r + " is " + missingNum)
	 */



	public static void NotNaive(int[][] s)
	{
		int sum = s.length * ((s.length +1)/2);

		for(int i = 0; i < s.length; i++)
		{
			int missing = sum;

			for(int j = 0; j < s.length; j++)
			{
				for(int p = 1; p <= s.length; p++)
				{
					if(s[i][j] == p);
					{
						missing -= s[i][j];
					}
				}
			}
			System.out.println("The missing number in row " + i + " is: " + missing);
		}
	}


	public static boolean headerprinted = false;


	public static void main(String[] args) throws Exception 
	{
		
		
		
		
		//GroupingTargets.createCSVTestData();
		//GroupingTargets.testWeka();

		int lstart=0,  lend=2, hstart=8, hend=10;

		int perc = 10;//Integer.parseInt(args[3]);
		/*int nrow=2;//Integer.parseInt(args[0]), 
		int ncol=2;//Integer.parseInt(args[1]);
*/		//int ITER=1;
		//int nTargets=nrow*ncol;
		//int nRes=1;
		//double dmax = 4;//Integer.parseInt(args[2]);

		int[][] e1 = {
				/*{30, 30, 85, 150},*/
				/*{40, 25, 125, 165},*/
				//{40, 30, 145, 200}
				/*{40, 40, 120, 266}*/
				/*50,40,185,333*/
		};
		
		
		int nrow = 14;
		
		int ncol = 14;
		
		int dmax = 56;
		
		int graphk = 33; // number of cluster when I built an example
		
		int ITER = 5;
		
		int als[] = {2}; //DO + weka + CON target per cluster
		//radius
		int rs[] = {1}; // PACMAN
		int rs1[] = {1}; // attack cluster
		int rs2[] = {1}; // attack cluster + weka
		int percentagethreshold = 30; // perc of targets in restricted set other than attack set
		int als1[] = {2};// #targets in cluster for Attackcluster + weka
		
		int ranges[][] = {{0,1},{3,8},{9, 10}};
		int[] percforranges = {75, 20, 5};
		
		/*int ranges[][] = {{0,7},{6,8},{8, 10}};
		int[] percforranges = {90, 0, 10};*/
		
		int nRes=2;
		
		int solverk = 20; // number of cluster for solver
		
		int abstractionlevel = 2; // For DO with COntraction and weka clus
		
		int RADIUS = 1;
		
		
		
		
		
		
		int blockdim = 2; // block = blockdim x blockdim
		
		// nrow has to be divisible by block
		int naivencluster = (nrow*ncol)/(blockdim*blockdim);
		
		

		
		
		/*int nrow =Integer.parseInt(args[0]);
		int ncol = Integer.parseInt(args[1]);
		int dmax = Integer.parseInt(args[2]);
		int k = Integer.parseInt(args[3]);
		int RADIUS = Integer.parseInt(args[4]);*/
		
		
		
				
		
		int utiliy_l=0;
		int utility_h=10;
		int nTargets = nrow*ncol;
		
		int ncat = 3;
		
		int[] targetsincat = getTargetsInCats(nTargets, percforranges);
		double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);
		
		
		int base = 0;
		int dest = 0;
		
		int radius = 3;
		
		
		int dlim = 5;
		
		
		int ap = 4; // should be </= than cluster size
		
		/*int dmaxsuper = 30;
		int dminsuper = 5;*/
		
		
		
		
		
		
		// best algo 
		/*
		SecurityGameContraction.targets.clear();
		
		double[][] density=SecurityGameContraction.generateRandomDensity( perc, ITER, lstart, lend,  hstart, hend, nTargets, false);
		//double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);



		int percentages[] ={100};//{Integer.parseInt(args[1])};
		int thresholds[] = {2};
		
		SecurityGameContraction.doubleOracleTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		*/
		
		
		
		///////////////////////
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		HashMap<Integer, ArrayList<TargetNode>> alltargets = new HashMap<Integer, ArrayList<TargetNode>>();
		HashMap<Integer, HashMap<Integer, TargetNode>> alltargetmaps = new HashMap<Integer, HashMap<Integer, TargetNode>>();
		//HashMap<Integer, ArrayList<Integer>[]> allclus = new HashMap<Integer, ArrayList<Integer>[]>();
		HashMap<Integer, ArrayList<Integer>[]> allclus = new HashMap<Integer, ArrayList<Integer>[]>();
		//double[][] density=SecurityGameContraction.generateRandomDensity( perc, ITER, lstart, lend,  hstart, hend, nTargets, false);
		
		//double[][] density = new double[ITER][nTargets];
		
		
		
		for(int iter = 0; iter<ITER; iter++)
		{
			ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
			HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
			ClusterTargets.buildcsvGraphExp(nrow,ncol,density,targets, iter );
			//SecurityGameContraction.assignRandomDensityZeroSum(density, gamedata, targets, iter);
			//SecurityGameContraction.buildGraph(nrow, ncol, gamedata, targets);
			//SecurityGameContraction.assignRandomDensityZeroSum(density, gamedata, targets, iter);
			alltargets.put(iter, targets);
			for(TargetNode t : targets)
			{
				targetmaps.put(t.getTargetid(), t);
				
			}
			alltargetmaps.put(iter, targetmaps);
			ArrayList<Integer>[] clus = GroupingTargets.makeClusterWithRange(graphk , nTargets, utiliy_l, utility_h, 
					targets, targetmaps, density, iter, ranges, percforranges);
			ClusterTargets.buildFile(nrow,ncol,density,targets, iter );
			allclus.put(iter, clus);
			
			
		}
		
		
		
		
		
		
		
		
		
		
		
		//////real world data/////
		
		
		
		int k=20;
		nrow = 50;
		ncol = 50;
		ITER = 1;
		GroupingTargets.wekaClusteringWithDOExpRW(nrow,ncol,base, dest, k, radius, dmax, nRes);
		
		SecurityGameContraction.DORWTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 10, 10 );
		
		
		int al = 4;
		
		ClusterTargets.dOWithWekaCONRWExp(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10, al);
		
		
		int r = 1;
		
		ClusterTargets.DOWithPACMANClusteringRWTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps,r , 10, 10);
		
		ClusterTargets.dOWithAttackClusterRWTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, r, 10, 10);
		
		
		int p = 30;
		
		ClusterTargets.dOWithAttackClusterAndWekaRWTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, r, 10, 10,  p, al);
		
		/////////////////////////////////////////
		
		
		
		
				//////////////exp1/////////////////
		
		
				//4 DO + GC multi + GP 3 + LP + GC multi 
				//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 10, 10 );
				SecurityGameContraction.targets.clear();
				
				for(int al: als)
				{
				
				// DO+ COntration + weka
				   //ClusterTargets.dOWithWekaCONExp(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10, al);
				}
				
				for(int r: rs)
				{
				
					//ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps,r , 10, 10);
				}
				
				for(int r: rs1)
				{
					//activate clustering from the beginning, slave limit 10, path 10
				   //ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, r, 10, 10);
				}
				
				//int all[] = {10,15,20,25,30,35,40,50};
				
				for(int r: rs2)
				{
				
					for(int al: als1)
					{
						
						int[] pt = {1,2,5,10,15,20,25,30};
						for(int p: pt)
						{
							//ClusterTargets.dOWithAttackClusterAndWekaTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, r, 10, 10,  p, al);
						}
					}
				}
		
		
				
				
				///////exp2 ///// effect of targets per cluster
				/*int al[] = {2,5,7,10,12};
				
				for(int a: al)
				{
				
				// DO+ COntration + weka
				   ClusterTargets.dOWithWekaCONExp(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10, a);
				}
				
				
				
				for(int a: al)
				{
					
					int[] pt = {30};
					for(int p: pt)
					{
						// r =2
						ClusterTargets.dOWithAttackClusterAndWekaTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 2, 10, 10,  p, a);
					}
				}*/
				
				//////////////////////////////////////
				
				
				///////exp3 ////// effect of radius
				
				/*
				int rss [] = {1,2,3,4};
				for(int r: rss)
				{

					ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps,r , 10, 10);
				}

				for(int r: rss)
				{
					//activate clustering from the beginning, slave limit 10, path 10
					ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, r, 10, 10);
				}

				for(int r: rss)
				{
					int[] pt = {30};
					for(int p: pt)
					{
						// al =2
						ClusterTargets.dOWithAttackClusterAndWekaTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, r, 10, 10,  p, 2);
					}
				}*/
				
				///////////////////////////////
				
				
				////exp4 ////use exp1
				
				
			
				
				
				
				
				
				
				
				
				
				
				
				
		
			////  activate clustering from the beginning, slave limit 10, path 10
		
			
			
			
			
			//ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10, 5);
		
		
		
		
		
		
		
		

		
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 10, 10 );
		SecurityGameContraction.targets.clear();
		
		// DO+ COntration + weka
		//ClusterTargets.dOWithWekaCONExp(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10, abstractionlevel);
		
		
		
		
		
		//ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		
	////  activate clustering from the beginning, slave limit 10, path 10
		//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);

		
		//intra cluster path tries to cover all the targets rather than shortest path (224->272)
		// DO + Incremental attack clustering
		
		//activate clustering in the first iteration , slave limit 10, path 10
		//ClusterTargets.dOWithAttackClusterTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		
		
	
		
		//activate when iter>limit but still adding path by slave, slave limit 10, path 10
		//ClusterTargets.dOWithAttackClusterTest2(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		
		
		
		
		
		
		
		
		
			//  activate clustering from the beginning, slave limit 10, path 10
		// problem of targets not having any neighbor due to dominated targets
		// considering dominated targets as single ST has huge runtime
		// for later improvement...
			//ClusterTargets.dOWithAttackClusterTest4(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
			
		
		/*int[] sks = {nTargets, nTargets/2, nTargets/5, nTargets/10};
		
		for(int sk: sks)
		{
		
		
		// SO + weka
			ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, sk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10);
		}
		
		
		
		int[] abs = {1,2,5,10};
		
		for(int alevel: abs)
		{

		 //DO + weka
		 ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10, alevel);
		}*/
		
		
		// try without contraction
		
		SecurityGameContraction.targets.clear();
	
		
		//ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		
		
		
		
		
		
		
		
		/////exp 9 cluster refined to coarse
		
		
/*
		int solverks[] = {144,72,36,24,12};

		SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 10, 10 );




		for(int slverk: solverks)
		{
			ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, slverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10);
			SecurityGameContraction.targets.clear();
		}



		ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10);
		SecurityGameContraction.targets.clear();




		ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);

		ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		SecurityGameContraction.targets.clear();
*/




		
		
		
		/////exp 8 increasing distance results in worse solution quality
		
		/*int dm[] = {35, 40, 45, 50,55,60};
		
		for(int dx: dm)
		{
			SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dx, nRes, alltargets, alltargetmaps, 10, 10 );
			SecurityGameContraction.targets.clear();
		}
		
		
		
		
		for(int dx: dm)
		{
			ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dx, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
			SecurityGameContraction.targets.clear();
		}
		
		
		for(int dx: dm)
		{
			ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, solverk, radius, dx, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10);
			SecurityGameContraction.targets.clear();
		}
		
		
		
		for(int dx: dm)
		{
			ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dx, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10, 5);
			
		
			SecurityGameContraction.targets.clear();
		}
		
		
		
		for(int dx: dm)
		{
			ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dx, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
			
			SecurityGameContraction.targets.clear();
		}
		
		
		
		for(int dx: dm)
		{
			ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dx, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
			SecurityGameContraction.targets.clear();
		}
		*/
		
		
		
		
		/*for(int i=0; i<nrow; i++)
		{
			for(int j=0; j<ncol; j++)
			{
				System.out.print(new DecimalFormat("#000").format(((i)*nrow)+j)   + " ");
			}
			System.out.println();
		}*/ 
		
		///exp 7///////////////
		
		
		
		/*int slavelimits[] = {20};
		int pathlimits[] = {20};

		for(int slavelimit: slavelimits)
		{
			for(int pathlimit: pathlimits)
			{
				SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, slavelimit, pathlimit );
				SecurityGameContraction.targets.clear();
				//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, slavelimit, pathlimit);
			}
		}
		
		
		int slavelimits1[] = {10};
		int pathlimits1[] = {10};

		for(int slavelimit: slavelimits1)
		{
			for(int pathlimit: pathlimits1)
			{
				//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, slavelimit, pathlimit );
				SecurityGameContraction.targets.clear();
				ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, slavelimit, pathlimit);
			}
		}*/



		/*int slavelimits2[] = {5,10,15,20};
		int pathlimits2[] = {10};

		for(int slavelimit: slavelimits2)
		{
			for(int pathlimit: pathlimits2)
			{
				SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, slavelimit, pathlimit );
				SecurityGameContraction.targets.clear();
				//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, slavelimit, pathlimit);
			}
		}

		
		int slavelimits4[] = {5,10,15,20};
		int pathlimits4[] = {10};

		for(int slavelimit: slavelimits4)
		{
			for(int pathlimit: pathlimits4)
			{
				//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, slavelimit, pathlimit );
				SecurityGameContraction.targets.clear();
				ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, slavelimit, pathlimit);
			}
		}



		int slavelimits3[] = {5,10,15,20};
		int pathlimits3[] = {5,10,15,20};

		for(int slavelimit: slavelimits3)
		{
			for(int pathlimit: pathlimits3)
			{
				if(slavelimit==pathlimit)
				{
					SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, slavelimit, pathlimit );
					SecurityGameContraction.targets.clear();
					//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, slavelimit, pathlimit);
				}
			}
		}
		
		int slavelimits5[] = {5,10,15,20};
		int pathlimits5[] = {5,10,15,20};

		for(int slavelimit: slavelimits5)
		{
			for(int pathlimit: pathlimits5)
			{
				if(slavelimit==pathlimit)
				{
					//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, slavelimit, pathlimit );
					SecurityGameContraction.targets.clear();
					ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, slavelimit, pathlimit);
				}
			}
		}



*/

		
		
		
		
		
		
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 10, 10 );
		SecurityGameContraction.targets.clear();

		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 10, 20 );
		SecurityGameContraction.targets.clear();
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 20, 10 );
		SecurityGameContraction.targets.clear();
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, 20, 20 );
		SecurityGameContraction.targets.clear();


		/*int dm[] = {60, 65, 70, 75, 80};
		
		for(int dx: dm)
		{
			SecurityGameContraction.DOTest(density,ITER,nrow, ncol, dx, nRes, alltargets, alltargetmaps, 20, 20 );
			SecurityGameContraction.targets.clear();
		}
		*/
		
		
		
		
		
		
		//intra cluster path tries to cover all the targets rather than shortest path (224->272)
		// DO + Incremental attack clustering
		
		//activate clustering in the first iteration , slave limit 10, path 10
		//ClusterTargets.dOWithAttackClusterTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 20);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 20);
		SecurityGameContraction.targets.clear();
		
		
		
		
	
		
		//activate when iter>limit but still adding path by slave, slave limit 10, path 10
		//ClusterTargets.dOWithAttackClusterTest2(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest2(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 20);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest2(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest2(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 20);
		SecurityGameContraction.targets.clear();
		
		
		
		
		
		//  activate clustering from the beginning, slave limit 10, path 10
		//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 20);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.dOWithAttackClusterTest3(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 20);
		SecurityGameContraction.targets.clear();
		
		
		
		
		
		// SO + weka
		//ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10);
		// SO + weka
		//ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 20);
		// SO + weka
		//ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 20, 10);
		// SO + weka
		//ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 20, 20);
		
		

		
		
		// adapt number of cluster based on targets to cluster

		// DO + weka
		//ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 10);
		// DO + weka
	//	ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 10, 20);
		// DO + weka
		//ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 20, 10);
		// DO + weka
		//ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, 20, 20);

		
		//GroupingTargets.groupingWithDOExp(base, dest, graphk, radius, dmax, nRes, nTargets, ITER, ap, null,  alltargets, alltargetmaps);
		
		
		
		
		
		
	//	ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		SecurityGameContraction.targets.clear();
	//	ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 20);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.DOWithPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 20);
		SecurityGameContraction.targets.clear();
		
		
		
	//	ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 10);
		SecurityGameContraction.targets.clear();
	//	ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 10, 20);
		SecurityGameContraction.targets.clear();
	//	ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 10);
		SecurityGameContraction.targets.clear();
		//ClusterTargets.DOWithSplitPACMANClusteringTest(density,ITER,nrow, ncol, dmax, nRes, alltargets, alltargetmaps, RADIUS, 20, 20);
		SecurityGameContraction.targets.clear();
		
		
		
		
		
		
		
		
		
		
		
		// SO + weka
		//ClusterTargets.wekaClusteringWithSOExp(nrow,ncol,base, dest, solverk, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps);
		
		
		
		//SO + naiveClustering
		//ClusterTargets.naiveClusteringWithSOExp(nrow,ncol,base, dest, naivencluster, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps, blockdim);
		
		
		//14 baseline
		//SecurityGameContraction.noContractionNoColumnGenerationTest(density, ITER, nrow, ncol, dmax, nRes, alltargets, alltargetmaps );
		SecurityGameContraction.targets.clear();
		
		//SecurityGameContraction.noContractionWithColumnGenerationTest(density, ITER, nrow, ncol, dmax, nRes, alltargets, alltargetmaps);
		SecurityGameContraction.targets.clear();
		
		//SecurityGameContraction.noContractionWithColumnGenerationHeuTest(density, ITER, nrow, ncol, dmax, nRes, alltargets, alltargetmaps);
		SecurityGameContraction.targets.clear();

		
		
		
		//ClusterTargets.naiveClusetringTest();
		
		
		//
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		//GroupingTargets.wekaClusteringWithDOExpRW(nrow,ncol,base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps);
		
		
		
		
		
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(alltargets, alltargetmaps, ITER, nTargets , dmax, nRes);
		SecurityGameContraction.targets.clear();
		
		
		
		
		
		
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		
		/*ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		ArrayList<Integer>[] clus = GroupingTargets.makeGraph(k, radius, dlim , nTargets, utiliy_l, utility_h, ap, targets, targetmaps);*/
		//GroupingTargets.createGraph2(targets, targetmaps);	
		//ArrayList<Integer>[] clus = GroupingTargets.makeDummyCluster(k);
		
		
		//only one path per super targets
		// but add path incrementally
		//GroupingTargets.groupTargetBaseline3Exp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
				
		
		
		
		//only one path per super targets
		// with modified LP
		//GroupingTargets.groupTargetBaseline2ExpNewMILP(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
		
		
		
		//only one path per super targets
		//GroupingTargets.groupTargetBaseline2Exp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
		
		
		
		//GroupingTargets.groupTargetTest();
		//GroupingTargets.groupTargetBaselineExp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
		
		SecurityGameContraction.targets.clear();
		
		//SecurityGameContraction.baselineForGroupTest(ITER, nTargets, dmax, nRes, alltargets, alltargetmaps);
 		//SecurityGameContraction.targets.clear();
		
		
		
		
		
		


		 targetsincat = getTargetsInCats(nTargets, percforranges);
		//double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);


		/*double[][] density=SecurityGameContraction.generateRandomDensity( perc, ITER, lstart, lend,  hstart, hend, nTargets, false);
		//double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);



		int percentages[] ={100};//{Integer.parseInt(args[1])};
		int thresholds[] = {2};
		
		SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();*/




		//grouping test
		//SecurityGameContraction.GroupingTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();








		/*** REAL WORLD DATA TEST
		 * 
		 * 
		 * 
		 */


		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiRealWorldDataTest(1,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		/**
		 * END
		 */



		/****
		 * 
		 * 
		 * Exact solver
		 */
		int[][] gamedata = new int[nTargets][4];

		//SecurityGameContraction.contractionWithDoubleOracle(gamedata, nTargets, nRes, density, dmax, 0, 5, 5);



		//Exact DO + GC multi + GP 3 + LP + GC multi 
			//SecurityGameContraction.doubleOracleGCMultiExactLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		//SecurityGameContraction.doubleOracleGCMultiGP3LP_OPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		//SecurityGameContraction.doubleOracleGCMultiGP3LP_lexicoOPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		//SecurityGameContraction.doubleOracleGCMultiGP3LP_TOPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		//SecurityGameContraction.doubleOracleGCMultiGP3LP_modifiedOPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		//SecurityGameContraction.targets.clear();







		/**
		 * experiments
		 */

		//1 DO + GC single + GP multi + LP + GC multi 
		//	SecurityGameContraction.doubleOracleGCSingleGPMultiLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//2 DO + GC multi + GP multi + LP + GC multi 
		//    SecurityGameContraction.doubleOracleGCMultiGPMultiLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//3 DO + GC single + GP 3 + LP + GC multi 
			//SecurityGameContraction.doubleOracleGCSingleGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		//5 DO + GC single + GP multi + LP + sample path 
		//	SecurityGameContraction.doubleOracleGCSingleGPMultiLPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//6 DO + GC multi + GP multi + LP + sample path
		//z	SecurityGameContraction.doubleOracleGCMultiGPMultiLPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//7 DO + GC single + GP 3 + LP + sample path 
		//	SecurityGameContraction.doubleOracleGCSingleGP3LPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//8 DO + GC multi + GP 3 + LP + sample path 
		//	SecurityGameContraction.doubleOracleGCMultiGP3LPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		// 9 single oracle + greedy cover + greedy path
		//	SecurityGameContraction.contractionWithSingleOracleGreedyTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//10 single oracle + 2 targets + all paths
		//	SecurityGameContraction.contractionWithSingleOracleTest1(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		// 11 basic + GC Single + ins + All paths + LP
		//	SecurityGameContraction.basicGCSingleInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//12 basic + GC Multi + ins + All paths + LP 
		//	SecurityGameContraction.basicGCMultiInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//13 basic + GC Multi + ins + GP 3 + LP 
		//	SecurityGameContraction.basicGCMultiInstantGP3LPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//14
		//SecurityGameContraction.noContractionNoColumnGenerationTest(density, ITER, nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//15
		//	SecurityGameContraction.noContractionWithColumnGenerationTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();








		/**
		 * end of experiments
		 */



		/**
		 * 
		 * GP3 vs GP3 with APSP
		 */

		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.doubleOracleGCMultiGP3APSPLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);







		/**
		 * APSP TESTS
		 */


		//SecurityGameContraction.doubleOracleAPSPGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);



		//12 basic + GC Multi + ins + All paths + LP 
		//SecurityGameContraction.basicGCMultiInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		//SecurityGameContraction.targets.clear();
		//SecurityGameContraction.basicAPSPGCMultiInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);

		/////

		/**
		 * APSP TESTS ENDS
		 */





		//SecurityGameContraction.testPathGenerationv2();
		//SecurityGameContraction.basicAbstraction();
		//1
		//SecurityGameContraction.noContractionNoColumnGenerationTest(density, ITER, nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//2
		//SecurityGameContraction.noContractionWithColumnGenerationTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		//3 basic + seq + single res cover +  sample path + MILP
		//SecurityGameContraction.basicSeqAbstractionAlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//4 basic + seq + single res cover +  sample path + LP
		//SecurityGameContraction.basicSeqAbstractionWithLPGreedyCover1AlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//5 basic + seq + mult res cover +  sample path + LP
		//SecurityGameContraction.basicSeqAbstractionWithLPGreedyCoverMultResAlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);



		SecurityGameContraction.targets.clear();
		//6 basic + ins + single res cover +  sample path + MILP 
		//SecurityGameContraction.basicInstantAbstractionAlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//7 basic + ins + single res cover +  sample path + LP 




		//SecurityGameContraction.instantContractionTest(density, ITER, nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.instantContractionTestAPSP(density, ITER, nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.spanningTreeTest(density,ITER);
		//SecurityGameContraction.testPathGenerationV3(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		//SecurityGameContraction.MPI();
		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.contractionWithSingleOracleGreedyTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.contractionWithSingleOracleTest1(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();






		//9 Double oracle +  mult res cover + greedy path 2(initial paths) + MILP + path sample 
		//SecurityGameContraction.contractionWithDoubelOracleTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//10 Double oracle +  mult res cover + greedy path 2(initial paths) + MILP + greedycover 2 
		//SecurityGameContraction.contractionWithDoubelOracleGreedyCoverTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		SecurityGameContraction.targets.clear();

		SecurityGameContraction.targets.clear();
		//13 Double oracle +  mult res cover + greedy path 3(initial paths) + LP + greedycover 2 
		//SecurityGameContraction.contractionWithDoubelOracleGreedyCover3LPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		//14 Double oracle +  mult res cover + greedy path 3(initial paths) + LP + path sample  
		//SecurityGameContraction.contractionWithDoubelOracleGreedyPath3TestLP(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);


		/**
		 * do tests with extreme pruning
		 */


		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.basicSeqAbstractionAlgorithmWithExtreamPrunningTest(density, ITER,nrow, ncol, percentages, dmax);
		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.basicInstantAbstractionAlgorithmWithExtrmPruningTest(density, ITER,nrow, ncol, percentages, dmax);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.instantContractionWithExtreamPruningTest(density, ITER, nrow, ncol, percentages, thresholds, dmax);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.testPathGenerationWithExtreamPruningV3(density, ITER, nrow, ncol, percentages, thresholds, dmax);




		/*int[] contractionsizes = {5 };
		ITER = 1;
		//int x = ITER-1;
		double origcoinval = 0;
		double[] vals = new double[ITER];

		int targetssize = 15;
		int row =5;
		int col = 3;           
		for(int contractionsize : contractionsizes)
		{
			System.out.println("\n contraction size "+ contractionsize + "...");
			double sumcoins = 0.0;
			long sumcontracttime = 0;
			long sumsolvetime = 0;
			double nduplicatetargets = 0;
			double valentedge = 0;
			double constrsize = 0;


			for(int iter =0; iter<(ITER); iter++)
			{

				int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
				//SecurityGameAbstraction.convertToZeroSum(gamedata);
				Date start = new Date();
				long l1 = start.getTime();

				ArrayList<TargetNode> dominatednodes =  SecurityGameContraction.testContraction(targetssize, contractionsize, gamedata, 6, row, col);
				//System.out.println("\n iter "+ iter+ ", #dominated nodes "+ dominatednodes.size());
				Date stop = new Date();
				long l2 = stop.getTime();
				System.out.println("\n iter "+ iter+ " contraction size: "+dominatednodes.size());
				constrsize+= dominatednodes.size();


				long diff = l2 - l1;
				sumcontracttime += diff;

				for(TargetNode x: dominatednodes)
				{
					System.out.print(x.getTargetid() + " ");

				}
				System.out.println();
         		SecurityGameContraction.removePathsToDominatedNodes(SecurityGameContraction.targets,dominatednodes);
				SecurityGameContraction.copyInOrigTargets(SecurityGameContraction.targets);
				SecurityGameContraction.removeDominatedTargets(SecurityGameContraction.targets,dominatednodes);

				SecurityGameContraction.printNodesWithNeighborsAndPath(dominatednodes, SecurityGameContraction.targets);

				SecurityGameContraction.addVirtualBase(0,targetssize, SecurityGameContraction.targets);
				//int q = SecurityGameContraction.calculateEdgesWithCoin();
				HashMap<Integer, Integer> nodewithcoins= SecurityGameContraction.calculateNodesWithCoin(SecurityGameContraction.targets);
				SecurityGameContraction.transformToDirectedGraph(SecurityGameContraction.targets);
				//SecurityGameContraction.printNodesWithNeighborsAndPath(dominatednodes);

				//SecurityGameContraction.printDuplicateTargets();


				//sumcoins += LpSolverIP.solve(SecurityGameContraction.targets, dominatednodes, 9);
				Date start1 = new Date();
				l1 = start1.getTime();
				double vvv= MIPSolver3.solve(SecurityGameContraction.duplicatetargets, SecurityGameContraction.targets, nodewithcoins, dominatednodes, 6);

				if(contractionsize==0)
				{
					vals[iter] = vvv;
				}
				if((contractionsize>0) && (vvv>vals[iter]))
				{
					System.out.println("\n");
				}
				sumcoins += vvv;
				Date stop1 = new Date();
				l2 = stop1.getTime();
				long diff1 = l2 - l1;
				sumsolvetime += diff1;
				//TargetNode goal = CoinCollection.coinCollection(SecurityGameContraction.targets, dominatednodes);


				long secondInMillis = 1000;
				long minuteInMillis = secondInMillis * 60;
				long hourInMillis = minuteInMillis * 60;
				long dayInMillis = hourInMillis * 24;


				long elapsedDays = diff / dayInMillis;
				diff = diff % dayInMillis;
				long elapsedHours = diff / hourInMillis;
				diff = diff % hourInMillis;
				long elapsedMinutes = diff / minuteInMillis;
				diff = diff % minuteInMillis;
				long elapsedSeconds = diff / secondInMillis;
				//sumcoins += goal.getCoinvalue();

				//System.out.println("\niter"+ iter+", coin collected "+ goal.getCoinvalue());
				//System.out.println("Running Time "+ elapsedDays+" days, "+ elapsedHours+" hrs, "+ elapsedMinutes+" mins, "+ elapsedSeconds+" sec, " + diff + "ms" );
				//double coins = CoinCollection.printPath(goal);
				//System.out.println("\niter"+ iter+", coin collected path value "+ coins);
				//SecurityGameContraction.unimportanttargets.clear();
				nduplicatetargets += SecurityGameContraction.duplicatetargets.size();
				valentedge += nodewithcoins.size();
				SecurityGameContraction.targets.clear();
				SecurityGameContraction.duplicatetargets.clear();
				SecurityGameContraction.origtargets.clear();

			}
			double avgcoin = sumcoins/ITER;
			double avgcontracttime = sumcontracttime/ITER;
			double avgsolvetime = sumsolvetime/ITER;
			double avgndup = nduplicatetargets/ITER;
			valentedge = valentedge/ITER;
			constrsize = constrsize/ITER;



			writeInFile(contractionsize ,avgcoin, avgcontracttime, avgsolvetime, avgndup, valentedge, constrsize);
		}*/

		/*//int numCluster = 10;
		double[] numResources = {20};
		//int numTarget = 200;
		int ITER = 100;
		int[] clusternumbers = {40,20,10,5};
		for(int  numCluster: clusternumbers )
		{
			//double sumsolmmr = 0;
			//double summaxintrvlmmr = 0;
			double sumsolpoly = 0;
			double summaxintrvl = 0;
			double sumsolperfectmmr = 0;
			//double summaxintrvlperfectmmr = 0;

			double sumavg = 0;
			//double summaxintrvlpoly2 = 0;
			//double sumsolperfectmmr2 = 0;
			//double summaxintrvlperfectmmr2 = 0;


			boolean isequalcluster = false;

			for(int iter=0; iter<ITER; iter++)
			{
				//int[][] gamedata = SecurityGameAbstraction.createdummyData(numTarget,1,1);
				int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
				SecurityGameAbstraction.convertToZeroSum(gamedata);
				List<Integer>[] clusteredtargets = KmeanClustering.getBestclusteredTargets(numCluster, gamedata, isequalcluster);
				int[] clustersizes = getClusterSizes(clusteredtargets);
				int maxinterval = KmeanClustering.getMaxInterval(clusteredtargets, gamedata);
				StrategyMapping strmap = new StrategyMapping(clusteredtargets, numCluster, gamedata);

				int[][][] abstractgame = SecurityGameAbstraction.makeAbstractIntervalSecurityGame(strmap);
				double[][] avgabstractgame = SecurityGameAbstraction.makeAVGPayoffAbstractSecurityGame(strmap);
				double origexp = SecurityGameAbstraction.originalGameTesing(numResources, gamedata);
				//strmap.printSecurityGameMapping();
				//SecurityGameAbstraction.makePositive(gamedata, -125);
				//SecurityGameAbstraction.addNoise(gamedata, 5, 10, true);
		 *//**
		 * try with interval abstract game
		 *//*
				//double[] resultmmr = SecurityGameAbstraction.testingBoundedRationalMMR(numResources, gamedata,numCluster);
				//sumsolmmr += resultmmr[0];
				//summaxintrvlmmr += resultmmr[1];
				double[] resultperfectmmr = SecurityGameAbstraction.testingPerfectRationalMMR(numResources, 
						gamedata, abstractgame, strmap, numCluster, isequalcluster, origexp, clustersizes);
				sumsolperfectmmr += resultperfectmmr[0];
				summaxintrvl += maxinterval;

				double[] resultpoly = SecurityGameAbstraction.wildlifeAbstraction(numResources,
						gamedata, abstractgame, strmap, numCluster, isequalcluster, origexp, clustersizes);
				sumsolpoly += resultpoly[0];
				//summaxintrvlpoly += maxinterval;


		  *//**
		  * try with avg payoff abstract game
		  *//*

				double[] resultavgabsgm = SecurityGameAbstraction.testingAvgAbstraction(numResources, 
						gamedata, avgabstractgame, strmap, numCluster, isequalcluster, origexp, clustersizes);
				sumavg += resultavgabsgm[0];
				//summaxintrvlperfectmmr2 += maxinterval;

				double[] resultpoly2 = SecurityGameAbstraction.wildlifeAbstraction(numResources,
						gamedata, avgabstractgame, strmap, numCluster, isequalcluster, origexp);
				sumsolpoly2 += resultpoly2[0];
				//summaxintrvlpoly2 += maxinterval;
			}
			sumsolperfectmmr /= ITER;
			sumsolpoly /= ITER;
			sumavg /= ITER;
			summaxintrvl /= ITER;
			System.out.println("Sol quality mmr "+ sumsolperfectmmr+ " maxinterval "
			+ summaxintrvl + " sumsolpoly "+ sumsolpoly );
			writeInFile(numCluster, sumsolperfectmmr, sumsolpoly, sumavg, summaxintrvl);
		}
		//SecurityGameAbstraction.testingAvgAbstraction(numResources, gamedata, numCluster);
		   */
	}











	private static int[] getTargetsInCats(int nTargets, int[] percforcats) {


		int x[] = new int[percforcats.length];


		int sum = 0;

		for(int i=0; i<percforcats.length; i++)
		{
			x[i] = (int)Math.floor(nTargets*(percforcats[i]/100.00));
			sum += x[i];
		}




		if(sum<nTargets)
		{

			int max = Integer.MIN_VALUE;
			int index= 0;
			for(int i=0; i<percforcats.length; i++)
			{
				if(max>percforcats[i])
				{
					max = percforcats[i];
					index = i;
				}
			}

			x[index] += (nTargets-sum);

		}



		return x;
	}











	private static void writeInFile(int contractionsize, double avgcoin,
			double avgcontracttime, double avgsolvetime, double avgndup, double valentedge, double constrsize) 
	{

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(contractionsize+ "," + avgcoin  
					+ "," + avgcontracttime +","+avgsolvetime+
					","+avgndup+","+valentedge+","+constrsize+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}

	private static int[] getClusterSizes(List<Integer>[] clusteredtargets) {

		int[] clustersizes = new int[clusteredtargets.length];
		for(int i=0; i<clusteredtargets.length; i++)
		{
			clustersizes[i] = clusteredtargets[i].size();
		}
		return clustersizes;

	}

	/*public static void writeInFile(int numCluster, double sumsolmmr,
			double sumsopoly, double sumavg, double summaxintrvl, double valentedge, double constrsize) {
		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(numCluster+ "," + sumsolmmr  + "," + sumsopoly +"," +
			sumavg + ","+summaxintrvl+","+constrsize+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	 */




}
