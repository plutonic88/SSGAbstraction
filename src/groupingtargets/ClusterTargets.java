package groupingtargets;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;

import cs.Interval.ILP.MIPSolver4;
import cs.Interval.contraction.Logger;
import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;
import cs.com.allpair.AllPairShortestPath;
import weka.clusterers.EM;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;

public class ClusterTargets {
	
	
	public static Random rand1 = new Random(100);
	public static final int INFINITY = 99999;
	
	
	
	public static void printNodesWithNeighborsAndPath(  HashMap<Integer, TargetNode> targets) 
	{
		for(TargetNode node : targets.values())
		{
			System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");
			Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{
				
				
				for(TargetNode neighbor: node.getNeighbors())
				{
					
					if(node.getTargetid()==224 && neighbor.getTargetid()==272)
					{
						System.out.println();
					}
					
					System.out.println("---Neighbor : "+ neighbor.getTargetid());
					//Logger.logit("---Neighbor : "+ neighbor.getTargetid()+"\n");
					/**
					 * print path
					 */
					ArrayList<TargetNode> path = node.getPath(neighbor);
					System.out.print("Path : "+ node.getTargetid()+ " --> ");
					//Logger.logit("Path : "+ node.getTargetid()+ " --> ");
					for(TargetNode pathnode : path)
					{
						System.out.print(pathnode.getTargetid()+" --> ");
						Logger.logit(pathnode.getTargetid()+" --> ");
					}

					System.out.print(neighbor.getTargetid()+ "\n");
					System.out.println("Distance : " + node.getDistance(neighbor));
					//Logger.logit(neighbor.getTargetid()+ "\n");
					//Logger.logit("Distance : " + node.getDistance(neighbor));

					System.out.print("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");
					//Logger.logit("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");

				}
			}
		}

	}
	
	
	public static void buildcsvGraphExp(int numRow, int numCol, double[][] u,  ArrayList<TargetNode> targets, int iter) throws Exception 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);
		
		
		
		/*try {
			
			
			 File f = new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata"+iter+".csv");
			 
			 if(f.exists())
			 {
				 f.delete();
				 f.createNewFile();
			 }
			
			
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata"+iter+".csv"),true));
			
			pw.append("Id,U,X,Y"+"\n");
			pw.close();
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}*/

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, 0);
			/*node.defenderreward = gamedata[target][0];
			node.defenderpenalty = gamedata[target][1];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = gamedata[target][3];*/

			targets.add(node);
			if(target==0)
			{
				//graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;
		
		
		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				TargetNode tmp = targets.get(targetid);
				
				
				tmp.setRowCol(row, col);
				
				
				/*tmp.setCoinvalue(u[iter][targetid]);
				tmp.defenderreward = 0;
				tmp.defenderpenalty = -u[iter][targetid];
				tmp.attackerreward = u[iter][targetid];
				tmp.attackerpenalty = 0;
				tmp.setAnimaldensity(u[iter][targetid]);*/
				
				/*try {
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata"+iter+".csv"),true));
					
					pw.append(tmp.getTargetid()+","+u[iter][targetid]+ ","+(row*1) + ","+(col*1)+ "\n");
					pw.close();
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
*/
				
				
				
				
				
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						targets.get(targetid).addNeighbor(targets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						targets.get(targetid).setPath(targets.get(neighborid), pathnodes);
						targets.get(targetid).setPathUtility(targets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						//Double distance = 1.0;//rand.nextDouble()*10+5;
						
						double d1 = 1;
						
						/*System.out.println("target "+ targetid + ", neirow "+ neighborrow
								+ ", neicol "+ neighborcol + ", ele "+e[neighborrow][neighborcol]+", d = "+ d1);
						
						*/
						targets.get(targetid).addDistance(targets.get(neighborid), Math.floor(d1));
						//targets.get(neighborid).addDistance(targets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}

	}
	
	
	public static void buildFile(int numRow, int numCol, double[][] u,  ArrayList<TargetNode> targets, int iter) throws Exception 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);
		
		
		
		try {
			
			new File("result").mkdirs();
			
			 File f = new File("result/realdata"+iter+".csv");
			 
			 if(f.exists())
			 {
				 f.delete();
				 f.createNewFile();
			 }
			
			
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result/realdata"+iter+".csv"),true));
			
			pw.append("Id,U,X,Y"+"\n");
			pw.close();
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;
		
		
		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				TargetNode tmp = targets.get(targetid);
				

				tmp.setCoinvalue(u[iter][targetid]);
				tmp.defenderreward = 0;
				tmp.defenderpenalty = -u[iter][targetid];
				tmp.attackerreward = u[iter][targetid];
				tmp.attackerpenalty = 0;
				tmp.setAnimaldensity(u[iter][targetid]);
				
				try {
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result/realdata"+iter+".csv"),true));
					
					pw.append(tmp.getTargetid()+","+u[iter][targetid]+ ","+(row*1) + ","+(col*1)+ "\n");
					pw.close();
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

				targetid++;
			}
		}

	}


	
	
	
	public static void wekaClusteringWithDOExp(int nrow, int ncol, int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int slavelimit, int pathlimit ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;
		long clusteringtime = 0;
		long slavetime = 0;
		int finalsize = 0;
		int totalslaveiter = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			//ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			//printNodesWithNeighborsAndPath(targetmaps);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = wekaClusteringWithDO(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, iter, slavelimit, pathlimit);
			//double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime};
			
			
			totalslaveiter += res[7];
			sumdefexp += res[0];
			solvingtime+= res[2];
			revmaptime += res[6];
			clusteringtime += res[1];
			slavetime += res[5];
			finalsize += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		finalsize /= LIMIT;
		clusteringtime /= LIMIT;
		slavetime /= LIMIT;
		totalslaveiter /= LIMIT;
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		//writeInFileST("DOWithWeka",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime, nTargets);
		
		SecurityGameContraction.writeInFile("DOWithWeka",(int)finalsize, sumdefexp, 0,solvingtime, slavetime,totaltime, nTargets, totalslaveiter,clusteringtime, slavelimit, pathlimit);
		
		

	}
	
	
	public static void wekaClusteringWithSOExp(int nrow, int ncol, int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int slavelimit, int pathlimit ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;
		long clusteringtime = 0;
		long slavetime = 0;
		int finalsize = 0;
		int totalslaveiter = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			//ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			//printNodesWithNeighborsAndPath(targetmaps);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = wekaClusteringWithSO(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, iter, nrow, ncol, slavelimit, pathlimit);
			//double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime};
			
			
			sumdefexp += res[0];
			solvingtime+= res[2];
			revmaptime += res[6];
			totalslaveiter += res[7];
			clusteringtime += res[1];
			slavetime += res[5];
			finalsize += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		finalsize /= LIMIT;
		clusteringtime /= LIMIT;
		slavetime /= LIMIT;
		totalslaveiter /= LIMIT;
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		//writeInFileST("DOWithWeka",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime, nTargets);
		
		SecurityGameContraction.writeInFile("SOWithWeka",(int)finalsize, sumdefexp, 0,solvingtime, slavetime,totaltime, nTargets, totalslaveiter,clusteringtime, slavelimit, pathlimit);
		
		

	}
	
	
	public static void naiveClusteringWithSOExp(int nrow, int ncol, int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int blocksize ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;
		long clusteringtime = 0;
		long slavetime = 0;
		int finalsize = 0;
		int totalslaveiter = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			//ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			//printNodesWithNeighborsAndPath(targetmaps);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = naiveClusteringWithSO(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, iter, nrow, ncol, blocksize);
			//double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime};
			
			
			sumdefexp += res[0];
			solvingtime+= res[2];
			revmaptime += res[6];
			totalslaveiter += res[7];
			clusteringtime += res[1];
			slavetime += res[5];
			finalsize += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		finalsize /= LIMIT;
		clusteringtime /= LIMIT;
		slavetime /= LIMIT;
		totalslaveiter /= LIMIT;
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		//writeInFileST("DOWithWeka",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime, nTargets);
		
		SecurityGameContraction.writeInFile("SOWithNaive",(int)finalsize, sumdefexp, 0,solvingtime, slavetime,totaltime, nTargets, totalslaveiter,clusteringtime);
		
		

	}
	
	
	
	

	private static void writeInFileST(String algo, int finalsize, double defexp, long solvingtime, long revmaptime, long clusteringtime, long totaltime, int nTargets) 
	{
		
		//ClusteringWithDO",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"grp-result.csv"),true));
			pw.append(algo+ ","+nTargets+","+finalsize+","+defexp+"," + clusteringtime+ ","+solvingtime+ ","+revmaptime+","+totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	
	
	
	public static void buildAPSPWeka(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int nRes,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int[][] apspmat, AllPairShortestPath apsp ) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		// map = new HashMap<Integer, Integer>();
		// mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		
		 int[][] apspmat1 =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		SecurityGameContraction.purifyAPSPMatrixZeroGT(apspmat1, targets, nTargets, map, mapback);
		
		//apspmat = new int[apspmat1.length][apspmat1[0].length];
		
		
		
		for(int i=0; i<apspmat1.length; i++)
		{
			for(int j=0; j<apspmat1[0].length; j++)
			{
				apspmat[i][j] = apspmat1[i][j];
			}
		}




		return;
	}
	
	
	public static ArrayList<Integer> buildGreedyCoverMultResWeka(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int nRes,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int[][] apspmat, AllPairShortestPath apsp ) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		// map = new HashMap<Integer, Integer>();
		// mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		
		 int[][] apspmat1 =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		SecurityGameContraction.purifyAPSPMatrixZeroGT(apspmat1, targets, nTargets, map, mapback);
		
		//apspmat = new int[apspmat1.length][apspmat1[0].length];
		
		
		
		for(int i=0; i<apspmat1.length; i++)
		{
			for(int j=0; j<apspmat1[0].length; j++)
			{
				apspmat[i][j] = apspmat1[i][j];
			}
		}



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);

		tcur =  SecurityGameContraction.greedyCoverMultRes(base, targets, dmax, targetssorted, apspmat1, map,mapback, nRes);

		return tcur;
	}
	
	
	private static void makeAdjacencyMatrix(int[][] adjacencymatrix,
			ArrayList<TargetNode> targetmaps, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		int i=1; 
		for(TargetNode n: targetmaps)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				//System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();
				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}





	}
	
	
	
	private static void makeAdjacencyMatrixST(int[][] adjacencymatrix,
			ArrayList<TargetNode> targetmaps, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback, HashMap<Integer,TargetNode> nodes) {


		int i=1; 
		for(TargetNode n: targetmaps)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				//System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");


				if(nodes.keySet().contains(nei.getTargetid()))
				{
					adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
				}
			}
			i++;
		}


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();
				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}





	}


	
	
	private static void printSuperTargets(HashMap<Integer, SuperTarget> sts, HashMap<Integer,ArrayList<Integer>> stpaths, HashMap<Integer,Double> dstravel) {



		Iterator<SuperTarget> itr = sts.values().iterator();

		//for(SuperTarget st : itr)
		while(itr.hasNext())
		{
			SuperTarget node = itr.next();
			System.out.println("\n\n******Super target node " + node.stid+"******");
			//Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				// print the nodes
				System.out.print("---Nodes : ");
				for(TargetNode n: node.nodes.values())
				{
					System.out.print( n.getTargetid()+ " ");

				}
				if(node.nodes.size()>1 /*&& stpaths.keySet().contains(node.stid)*/)
				{
					System.out.print("\n---Travelpath : ");
					for(Integer n: stpaths.get(node.stid))
					{
						System.out.print(n+"->");
					}
					System.out.print("\n---TravelDist : " + dstravel.get(node.stid));
					
				}
				System.out.print("\n---Neighbors : ");
				for(SuperTarget neighbor: node.neighbors.values())
				{
					System.out.print(neighbor.stid+ "("+node.distances.get(neighbor)+"), ");

				}
				ArrayList<TargetNode> already = new ArrayList<TargetNode>();
				System.out.print("\n---Neighbor Targets : ");
				for(TargetNode neighbor: node.ap.values())
				{
					for(TargetNode nei: neighbor.getNeighbors())
					{
						if(!already.contains(nei) && !node.nodes.values().contains(nei))
						{
							System.out.print(nei.getTargetid()+ " ");
							already.add(nei);
						}
					}

				}
				System.out.print("\n---AP : ");
				for(TargetNode a: node.ap.values())
				{
					System.out.print(a.getTargetid()+ " ");

				}
			}
		}

	}
	
	
	private static void printSuperTargets(HashMap<Integer, SuperTarget> sts) {



		Iterator<SuperTarget> itr = sts.values().iterator();

		//for(SuperTarget st : itr)
		while(itr.hasNext())
		{
			SuperTarget node = itr.next();
			System.out.println("\n\n******Super target node " + node.stid+"******");
			//Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				// print the nodes
				System.out.print("---Nodes : ");
				for(TargetNode n: node.nodes.values())
				{
					System.out.print( n.getTargetid()+ " ");

				}
				System.out.print("\n---Neighbors : ");
				for(SuperTarget neighbor: node.neighbors.values())
				{
					System.out.print(neighbor.stid+ " ");

				}
				ArrayList<TargetNode> already = new ArrayList<TargetNode>();
				System.out.print("\n---Neighbor Targets : ");
				for(TargetNode neighbor: node.ap.values())
				{
					for(TargetNode nei: neighbor.getNeighbors())
					{
						if(!already.contains(nei) && !node.nodes.values().contains(nei))
						{
							System.out.print(nei.getTargetid()+ " ");
							already.add(nei);
						}
					}

				}
				System.out.print("\n---AP : ");
				for(TargetNode a: node.ap.values())
				{
					System.out.print(a.getTargetid()+ " ");

				}
			}
		}

	}

	
	
	
	private static void assignSTValues(HashMap<Integer, SuperTarget> currentst,
			HashMap<Integer, TargetNode> targetmaps) {

		
		
		for(SuperTarget st: currentst.values())
		{
			double maxval = Double.MIN_VALUE;
			for(TargetNode t: st.nodes.values())
			{
				if(t.attackerreward>maxval)
				{
					maxval = t.attackerreward;
					st.attackerreward = t.attackerreward;
					st.attackerpenalty = 0;
					
					st.defenderreward = 0;
					st.defenderpenalty = -t.attackerreward;
				}
			}
		}
		
	}



	static HashSet combine(Integer[] arr, int k, int startId, int[] branch, int numElem,HashSet arrSet)
	{
		if (numElem == k)
		{
			//System.out.println("k: "+k+(Arrays.toString(branch)));
			ArrayList<Integer> mySet = new ArrayList<Integer>();
			for(int i=0;i<branch.length;i++)
			{
				mySet.add(branch[i]);
			}
			arrSet.add(mySet);
			return arrSet;
		}

		for (int i = startId; i < arr.length; ++i)
		{
			branch[numElem++]=arr[i];
			combine(arr, k, ++startId, branch, numElem, arrSet);
			--numElem;
		}
		return arrSet;
	}

	
	
	private static int findIndex(Instances newinstance, Integer tid) {
		// TODO Auto-generated method stub
		
		int index = 0;
		
		Iterator<Instance> ins = newinstance.iterator();
		
		while(ins.hasNext())
		{
			Instance newins = ins.next();
			if(newins.value(0) == tid.intValue())
				return index;
			index++;
		}
		
		
		return -1;
	}
	
	

	private static boolean isNeighbor(SuperTarget st1, SuperTarget st2) {
		
		
		// if they have nodes in neighbor
		if(st1.neighbors.values().contains(st2))
			return true;
		
		return false;
	}
	
	
	
private static boolean areBothNei(SuperTarget s1, SuperTarget s2, SuperTarget tempst) {
		
		
		if(isNeighbor(s1, tempst) && isNeighbor(s2, tempst))
			return true;
		
		return false;
	}

	
	public static HashMap<Integer, SuperTarget> clusterTargetsWeka(ArrayList<Integer> targetstocluster, 
			ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps, double dmax, int k, int radius,
			HashMap<Integer, Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, EM dc,
			Instances instances, HashMap<Integer,Integer> apspmap, int[][] apspmat, AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback) throws Exception
	{
		
		
		
		
		
		/**
		 * make new instances with targets to cluster
		 */
		
		Instances newinstance = new Instances(instances);
		
		
		
		
		
		/**
		 * remove unwanted instance
		 * index is the target id
		 */
		
		for(Integer tid: targetmaps.keySet())
		{
			if(!targetstocluster.contains(tid))
			{
				
				
				
				
				int index = findIndex(newinstance, tid);
				
				if(index != -1)
				{
					newinstance.remove(index);
				}
				
			}
		}
		
		
		System.out.println("New instance size  "+ newinstance.size());
		
		/**
		 * cluster
		 */
		
		int totalcluster = k + (instances.size()-newinstance.size());
		
		if(!targetstocluster.contains(0))
		{
			System.out.println("No base mmmmmmm  "+ newinstance.size());

		}
		
		ArrayList<Integer>[] clusters = clusterWithWeka(k, newinstance, totalcluster, targetstocluster,targetmaps, dc, apspmap, apspmat);
		
		
		//ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[targetmaps.size()];
		
		/*for(int t=0; t<targetmaps.size(); t++)
		{
			clusters[t] = new ArrayList<Integer>();
			clusters[t].add(t);
		}
		*/
		

		ArrayList<Integer> emptyindex = new ArrayList<Integer>();
		
		
		for(int i=0; i<clusters.length; i++)
		{
			if(clusters[i].size()==0)
			{
				totalcluster--;
				emptyindex.add(i);
			}
			
		}
		
		
		
		ArrayList<Integer>[] realclusters = (ArrayList<Integer>[])new ArrayList[totalcluster];
		
		for(int i=0; i<totalcluster; i++)
		{
			realclusters[i] = new ArrayList<Integer>();
		}
		
		
		int j=0;
		for(int i=0; i<clusters.length; i++)
		{
			if(clusters[i].size()!=0)
			{
				for(Integer n: clusters[i])
				{
					realclusters[j].add(n);
				}
				j++;
			}
			
		}
		
		
			
		
		
		
		
		
		HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(realclusters,targetmaps);
		
		/**
		 * for every super target determine 2 AP
		 */
		
		
		
		//printClusters(realclusters);
		



		// printSuperTargets(sts);
		int stssize = 1; // to keep track if any targets were clustered
		
		/**
		 * 2. For every pair of cluster: minimize the distance d : d1 + d(a1,a2) + d2 
		 * 3. Find a1 and a2 which will minimize the distance
		 */

		/**
		 *  Make a temporary merged supertarget
		 *  Then for every pair of access points
		 *  try to find the minimum one
		 *  save it the config for a particular d
		 *  
		 */

		System.out.println("Computing dmin");
		
		
		//chooseAP(sts, apsp, apspmat, apspmap, apspmapback, dstravel, targetmaps, dmax, stpaths);
		
		
		//printNodesWithNeighborsAndPath(targetmaps);
		
		for(int i = 0; i<sts.size(); i++)
		{
			SuperTarget st = sts.get(i);
			
			

			// if cluster is chnaged then compute AP
			// if not changed then use the previous AP
			
			chooseDegreeBasedAP(sts, targetmaps, sts.get(i), apsp, apspmat, apspmap, apspmapback, dstravel, stpaths, (int)dmax);
			
			
		}
		
		
		

		//}
		//updateNeighbors(sts);
		
		//printNodesWithNeighborsAndPath(targetmaps);
		//printSuperTargets(sts);
		
		 return sts;
	}
	
	
	public static HashMap<Integer, SuperTarget> clusterTargetsNaive(ArrayList<Integer> targetstocluster, 
			ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps, double dmax, int k, int radius,
			HashMap<Integer, Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths,
			HashMap<Integer,Integer> apspmap, int[][] apspmat, AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback, int row, int blocksize) throws Exception
	{
		
		
		
		
		
		/**
		 * make new instances with targets to cluster
		 */
		
		
		
		ArrayList<Integer>[] clusters = naiveCluster(row,blocksize,k);
		printClusters(clusters);
		
		
		
		HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,targetmaps);
		
		/**
		 * for every super target determine 2 AP
		 */
		
		
		
		
		



		 printSuperTargets(sts);
		int stssize = 1; // to keep track if any targets were clustered
		
		/**
		 * 2. For every pair of cluster: minimize the distance d : d1 + d(a1,a2) + d2 
		 * 3. Find a1 and a2 which will minimize the distance
		 */

		/**
		 *  Make a temporary merged supertarget
		 *  Then for every pair of access points
		 *  try to find the minimum one
		 *  save it the config for a particular d
		 *  
		 */

		System.out.println("Computing dmin");
		
		
		//chooseAP(sts, apsp, apspmat, apspmap, apspmapback, dstravel, targetmaps, dmax, stpaths);
		
		
		
		
		for(int i = 0; i<sts.size(); i++)
		{
			SuperTarget st = sts.get(i);
			
			

			// if cluster is chnaged then compute AP
			// if not changed then use the previous AP
			
			
			
			
				
			
				chooseDegreeBasedAP(sts, targetmaps, sts.get(i), apsp, apspmat, apspmap, apspmapback, dstravel, stpaths, (int)dmax);
			
			
		}
		
		
		

		//}
		//updateNeighbors(sts);
		
		//printNodesWithNeighborsAndPath(targetmaps);
		//printSuperTargets(sts);
		
		 return sts;
	}
	

	
	

	
	private static ArrayList<Integer>[] naiveCluster(int nrow, int blockdim, int ncluster) {
		
		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[ncluster];
		
		for(int i=0; i<ncluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}
		
		
		for(int c=0; c<ncluster; c++)
		{
			
			
			int starttarget = 0;
			
			if(c==0) 
			{
				starttarget = 0;
			}
			else if(c%(nrow/blockdim)==0)
			{
				starttarget = c*blockdim*blockdim;
			}
			else
			{
				// get the k-1th targetid for previous cluster c-1
				int id = clusters[c-1].get(blockdim-1);
				starttarget = id+1;
			}
			
			
			int row = starttarget/nrow;
			int col = starttarget % nrow;
			
			
			for(int i=row; i<(row+blockdim); i++)
			{
				for(int j=col; j<(col+blockdim); j++)
				{
					int id = i*nrow + j;
					clusters[c].add(id);
				}
			}
			
			
			
		}
		
		return clusters;
		
	}


	private static void chooseAP(HashMap<Integer, SuperTarget> sts, AllPairShortestPath apsp,
			int[][] apspmat, HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback, 
			HashMap<Integer,Double> dstravel, HashMap<Integer,TargetNode> targetmaps, double dmax, HashMap<Integer,ArrayList<Integer>> stpaths) {
		

		for(SuperTarget tempst: sts.values())
		{
			//SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
			//tempst.stid = 200 + st1.stid + st2.stid;
			// printSuperTarget(tempst);
			/**
			 * for every pair of access points
			 */
			
			System.out.println("Computing ap for ST "+ tempst.stid);
			
			
			double mindi = Double.MAX_VALUE;
			int stid1 = -1;
			//int stid2 = -1;
			int aid1 = -1;
			int aid2 = -1;
			int s1id = -1;
			int s2id = -2;
			int s1ap = -1;
			int s2ap = -1;
			double sda1a2 = -1;
			double sd1 = -1;
			double sd2 = -1;
			ArrayList<Integer> a1a2path = new ArrayList<Integer>();
			ArrayList<Integer> a2path = new ArrayList<Integer>();
			ArrayList<Integer> a1path = new ArrayList<Integer>();
			SuperTarget dminst = new SuperTarget();
			
			double mindist = 500000;
			
			
			
			

			if(tempst.nodes.size()>1)
			{
				
				if(tempst.stid==14)
				{
					System.out.println("shortestdist(a1,s1) ");
				}



				//check if 0 is its neighbor if so just keep the ap that connects to base

				if(tempst.neighbors.containsKey(0))
				{
					//compute every parameter

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if(a1.getTargetid() != a2.getTargetid())
							{
								// for every pair of supertargets which are not st1 or st2
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if(s1.stid != s2.stid)
									{


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2


											double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);





											//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

											if(dista1a2 == 0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

												dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
											}

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}

				else
				{
					
					System.out.println("Here for ST" + tempst.stid);

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if(a1.getTargetid() != a2.getTargetid())
							{
								// for every pair of supertargets which are not st1 or st2
								for(SuperTarget s1: tempst.neighbors.values()) // need neighbor cluster of a1
								{
									for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
									{

										//if(s1.stid != s2.stid)
										{


											boolean arebothnei = areBothNei(s1,s2,tempst);

											if(s1.stid != tempst.stid && 
													s2.stid != tempst.stid  && arebothnei)
											{
												// measure the distance between a1->s1 and a2->s2

												//a1->s1 distance between a target and a supertarget
												/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


												ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
												ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
												ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


												double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
												//TODO
												//System.out.println("shortestdist(a1,s1) "+ d1);

												double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
												//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
												//TODO
												//System.out.println("shortestdist(a2,s2) "+ d2);

												// next measure the intra cluster shortest traveling path using a1 and a2


												double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);





												//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

												if(dista1a2 == 0)
												{
													//throw new Exception("No path found to compute AP for st "+ tempst.stid);
													//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

													dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
												}

												//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

												// if any of the dist is <0 we know that it's not possible to have a path

												if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
												{
													double totaldi = d1[0] + dista1a2 + d2[0];
													if(totaldi < mindi)
													{
														//stid1 = .stid;
														//stid2 = st2.stid;
														dminst = tempst;
														aid1 = a1.getTargetid();
														aid2 = a2.getTargetid();
														mindi = totaldi;
														sda1a2 = dista1a2;
														s1id = s1.stid;
														s2id = s2.stid;
														sd1 = d1[0];
														sd2= d2[0];
														s1ap = (int)d1[1];
														s2ap = (int)d2[1];
														//spath.add(tmpspath.get(0));
														a1a2path.clear();
														for(Integer in: tmpa1a2spath)
														{
															a1a2path.add(in);
														}


														a1path.clear();
														for(Integer in: tmpa1path)
														{
															a1path.add(in);
														}


														a2path.clear();
														for(Integer in: tmpa2path)
														{
															a2path.add(in);
														}

														/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


													}
												}

											}
										}
									}
								}
							}
						}
					}
				}
			}

			//}
			//						}
			//
			//					}
			//				}
			//			}


			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
			//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			/**
			 * Merge two cluster which have min d
			 * For the new id of supertargets concat the strings with comma
			 * remove the old two clusters
			 * add the new one. 
			 */
			
			
			//TODO what to do when there is no path from a1 to a2?????
			
			
			
			if(sda1a2==-1 && tempst.nodes.size()>1)
			{
				System.out.println("shortestdist(a1,s1) ");
				//printSuperTarget(tempst);
			}

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				
				// Issue: How to fix the issue : connection to base node might be cut
				// Fix: Always keep the ap which connects to base node
				
				
				
				System.out.println("AP done for st "+ tempst.stid);
				
				dstravel.put(tempst.stid, sda1a2);
				stpaths.put(tempst.stid, a1a2path);
				//SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
				// printSuperTarget(newst);
				/*sts.remove(stid1);
				sts.remove(stid2);*/
				//printSuperTargets(sts);
				
				updateAP(tempst, sts, aid1, aid2);
				//sts.put(newst.stid, newst);
				//update the neighbors of ST
				// System.out.println("\n\n After merging # supertargets : "+ sts.size());
				// System.out.println("\n After merging new supertargets : ");
				
				/**
				 * add s1 and s2 as neighbor of st if they are not already
				 */
				
				addNei(tempst.stid, aid1, s1id, sd1, a1path, sts, targetmaps, s1ap);
				addNei(tempst.stid, aid2, s2id, sd2, a2path, sts, targetmaps, s2ap);
				addNei(aid1,aid2,sda1a2, a1a2path, tempst, targetmaps);
				
				
				/**
				 * add a1 a2 as neighbor if they are not already
				 */
				
				
				//printSuperTargets(sts);
				//printNodesWithNeighborsAndPath(targetmaps);
				
				
				//System.out.println("hi");
				 
			}
			
		}

		
		
	}


	private static void addNei(int aid1, int aid2, double sda1a2, ArrayList<Integer> a1a2path, SuperTarget tempst, HashMap<Integer,TargetNode> targetmaps) {
		
		
		//printSuperTarget(tempst);
		
		
		TargetNode srcnode = tempst.nodes.get(aid1);
		TargetNode destnode = tempst.nodes.get(aid2);
		//3, 11, 4, 13, 5, 12, 21, 20, 19
		
		if(!srcnode.getNeighbors().contains(destnode))
		{
			//TargetNode srcnode = sts.get(stid).nodes.get(aid1);
			//TargetNode destnode = sts.get(s1id).nodes.get(a1path.get(a1path.size()-1));
			
			if(!srcnode.getNeighbors().contains(destnode))
			{
				srcnode.addNeighbor(destnode);
				destnode.addNeighbor(srcnode);
				
				// add path
				
				
				ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
				
				for(Integer n: a1a2path)
				{
					if(n != srcnode.getTargetid() && n!= destnode.getTargetid())
					{
						pathnodes.add(targetmaps.get(n));
					}
				}
				
				srcnode.setPath(destnode, pathnodes);
				srcnode.addDistance(destnode, sda1a2);
				
				
				ArrayList<TargetNode> revpathnodes = new ArrayList<TargetNode>();
				
				for(int i=0; i<pathnodes.size(); i++)
				{
					revpathnodes.add(pathnodes.get(pathnodes.size()-i-1));
				}
				
				destnode.setPath(srcnode, revpathnodes);
				destnode.addDistance(srcnode, sda1a2);
				
				
			}
		}
		
		
		
		
	}


	private static void addNei(int stid, int aid1, int s1id, double sd1, ArrayList<Integer> a1path,
			HashMap<Integer, SuperTarget> sts, HashMap<Integer,TargetNode> targetmaps, int s1ap) {
		
		
		/**
		 * a1......s1....
		 */
		
		// first check if s1 is neibor
		// first check if s1 is neibor
		
		SuperTarget src = sts.get(stid);
		SuperTarget dest = sts.get(s1id);
		
		
		
		if(!src.neighbors.keySet().contains(dest.stid))
		{
			src.neighbors.put(dest.stid,dest);
			dest.neighbors.put(src.stid, src);
			
			// add path
			
			src.path.put(dest, a1path);
			src.distances.put(dest, sd1);
			
			
			ArrayList<Integer> revpath = new ArrayList<Integer>();
			
			for(int i=0; i<a1path.size(); i++)
			{
				revpath.add(a1path.get(a1path.size()-i-1));
			}
			
			dest.path.put(src, revpath);
			dest.distances.put(src, sd1);
			
			
			
			
		}
		

		
		
		// do the same thing for TargetNode
		
		
		
		TargetNode srcnode = sts.get(stid).nodes.get(aid1);
		TargetNode destnode = targetmaps.get(s1ap);
		
		if(!srcnode.getNeighbors().contains(destnode))
		{
			srcnode.addNeighbor(destnode);
			destnode.addNeighbor(srcnode);
			
			// add path
			
			
			ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
			
			for(Integer n: a1path)
			{
				if(n != srcnode.getTargetid() && n!= destnode.getTargetid())
				{
					pathnodes.add(targetmaps.get(n));
				}
			}
			
			srcnode.setPath(destnode, pathnodes);
			srcnode.addDistance(destnode, sd1);
			
			
			ArrayList<TargetNode> revpathnodes = new ArrayList<TargetNode>();
			
			for(int i=0; i<a1path.size(); i++)
			{
				revpathnodes.add(pathnodes.get(a1path.size()-i-1));
			}
			
			destnode.setPath(srcnode, revpathnodes);
			destnode.addDistance(srcnode, sd1);
			
			
			
			
		}
		
		
		
		
	}

	
	
private static void updateAP(SuperTarget curst, HashMap<Integer, SuperTarget> sts, int aid1, int aid2 ) {
		
		// remove old sts as neighbors
		
		
		
		//curst.neighbors.clear();
		curst.ap.clear();
		
		curst.ap.put(aid1, curst.nodes.get(aid1));
		curst.ap.put(aid2, curst.nodes.get(aid2));
		
		
		//update new neighbor for every st and curst
		for(SuperTarget st: sts.values())
		{
			if(curst.stid != st.stid)
			{
				if(!isPotentialNeighbor(curst, st))
				{
					st.neighbors.remove(curst.stid);
					st.distances.remove(curst);
					st.path.remove(curst);
					
					
					curst.neighbors.remove(st.stid);
					curst.distances.remove(st);
					curst.path.remove(st);
					
				}
			}
		}
		
	}


private static boolean isPotentialNeighbor(SuperTarget newst, SuperTarget st) {
	
	
	if(newst.stid == st.stid)
		return false;
	
	
	for(TargetNode nei: newst.ap.values())
	{
		for(TargetNode nei2 : st.ap.values())
		{
			if(nei.getNeighbors().contains(nei2))
				return true;
		}
	}
	
	return false;
}


	private static void printSuperTarget(SuperTarget sts) {



		//Iterator<SuperTarget> itr = sts.values().iterator();

		//for(SuperTarget st : itr)
		//while(itr.hasNext())
		//{
			SuperTarget node = sts;
			System.out.println("\n\n******Super target node " + node.stid+"******");
			//Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				// print the nodes
				System.out.print("---Nodes : ");
				for(TargetNode n: node.nodes.values())
				{
					System.out.print( n.getTargetid()+ " ");

				}
				System.out.print("\n---Neighbors : ");
				for(SuperTarget neighbor: node.neighbors.values())
				{
					System.out.print(neighbor.stid+ " ");

				}
				System.out.print("\n---AP : ");
				for(TargetNode a: node.ap.values())
				{
					System.out.print(a.getTargetid()+ " ");

				}
			}
		//}

	}

	
	
	
	
	private static double shortestdist(TargetNode a1, TargetNode a2, 
			HashMap<Integer,Integer> apspmap, int[][] apspmat, ArrayList<Integer> tmpa1path, 
			AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback, SuperTarget tempst, int dmax) {
		// TODO Auto-generated method stub



		double dmin = 0;
		//boolean isnei = false;
		ArrayList<Integer> minpath = new ArrayList<Integer>();
		
		//for(TargetNode t: s1.ap.values())
		//{


			//System.out.println("FInding shortest dist for target "+ a2.getTargetid());
			ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
			ArrayList<SuperTarget>  pnodes = new ArrayList<SuperTarget>();
			double distcovered = -1;
			if(a1.getNeighbors().contains(a2))
			{
				//pnodes = base.getPath(dest);
				int src = a1.getTargetid();
				int des = a2.getTargetid();

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];


			}
			else
			{




				int src = a1.getTargetid();
				int des = a2.getTargetid();


				

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];
				//System.out.print("dist covered "+ distcovered+"\n");



				ArrayList<Integer>	tmppathnodes = apsp.getPathInST(src, des, apspmap, apspmapback, tempst);

				for(int k=0; k<tmppathnodes.size(); k++)
				{
					pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
				}

				

				//throw new Exception("Base to not neighbor for initial set of paths **********8");

			}
			
			
			if(dmin<distcovered && (distcovered<=dmax))
			{

				dmin = distcovered;
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				//tmppath.add(a1.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
				//	System.out.print(dest.getTargetid()+"\n");
				//tmppath.add(a2.getTargetid());
				//System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				/*System.out.print("complete path : \n");
				for(int k=0; k<tmppath.size(); k++)
				{

					System.out.print(tmppath.get(k)+"->");
				}

				System.out.print("\n");
				*/
				minpath.clear();
				for(Integer p: tmppath)
				{
					minpath.add(p);
				}
			}

		//}
		
		
		
		tmpa1path.clear();
		for(Integer p: minpath)
		{
			tmpa1path.add(p);
		}
		
		
		
		return dmin;
	}
	
	
private static double shortestdist(TargetNode a1, TargetNode a2, SuperTarget tempst, double dmax, ArrayList<Integer> spath) {
		
		
	
	
	
	
	
		
		int nTargets = tempst.nodes.size(); 
		TargetNode start = new TargetNode(a1);
		//Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		ArrayList<TargetNode> closed = new ArrayList<TargetNode>();
		ArrayList<Integer> donetargets = new ArrayList<Integer>();
		
		
		int pathcounter = 0;
		int nodestocover = tempst.nodes.size();
		
		
		PriorityQueue<TargetNode> fringequeue = new PriorityQueue<TargetNode>(1000, new Comparator<TargetNode>() {  
		    
            public int compare(TargetNode w1, TargetNode w2) {                         
                return (int)(w1.distancecoveredyet - w2.distancecoveredyet);
            }      
        }); 
		
		fringequeue.add(start);
		
		
		
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==a2.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//System.out.println();
				//SecurityGameContraction.printPath(node);
				
				
				
				pathcounter++;
				//System.out.println();
				goals.clear();
				goals.add(node);
				spath.clear();
				SecurityGameContraction.makeClusterPathSeq(goals, spath);
				
				int count = coveredCount(spath);
				
				if(count == tempst.nodes.size())
				{
					return node.distancecoveredyet;
				}
				//break;

				/*if(pathcounter>5000)
					break;*/
			}
			else if(node.distancecoveredyet<=dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				//if(!isInclosed(closed,node))
				{
					closed.add(node);
					ArrayList<TargetNode> succs = ExpandTarget(node, tempst.nodes, dmax);
					/**
					 * add nodes to queue
					 */
					for(TargetNode suc: succs)
					{
						//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
						//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
						//if(!isInclosed(closed,node))
						{
							fringequeue.add(suc);
						}

					}
				}
				
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		
		
		return 0;
	}
	
	
	
	private static int coveredCount(ArrayList<Integer> spath) {
	
		int count = 0;
		
		ArrayList<Integer> done = new ArrayList<Integer>();
		
		for(Integer n: spath)
		{
			if(!done.contains(n))
			{
				done.add(n);
				count++;
			}
		}
		
		
		
	return count;
}


	public static boolean isInclosed(ArrayList<TargetNode> closed, TargetNode node) {

		
		for(TargetNode t: closed)
		{
			if(t.getTargetid() == node.getTargetid())
				return true;
		}
		return false;
		
	}


	private static ArrayList<TargetNode> ExpandTarget(TargetNode node, HashMap<Integer,TargetNode> nodes, double dmax ) 
	{
		ArrayList<TargetNode> successors = new ArrayList<TargetNode>();

		/**
		 * find the index for  node
		 */
		TargetNode tmpnode = nodes.get(node.getTargetid());
		
		for(TargetNode nei: tmpnode.getNeighbors())
		{
			if(nodes.values().contains(nei))
			{
				TargetNode newnei = new TargetNode(nei);
				newnei.distancecoveredyet = node.distancecoveredyet + tmpnode.getDistance(nei);
				
				newnei.parent = node;
				if(newnei.distancecoveredyet<=dmax)
				{
					successors.add(newnei);
				}
			}
		}
		return successors;

	}


	

	
	
	private static double[] shortestdist(TargetNode a1, SuperTarget s1, 
			HashMap<Integer,Integer> apspmap, int[][] apspmat, ArrayList<Integer> tmpa1path, 
			AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback) {
		// TODO Auto-generated method stub



		double dmin = Double.MAX_VALUE;
		//boolean isnei = false;
		ArrayList<Integer> minpath = new ArrayList<Integer>();
		int tid = -1;
		
		for(TargetNode t: s1.ap.values())
		{


			//System.out.println("FInding shortest dist for target "+ s1.stid);
			ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
			ArrayList<SuperTarget>  pnodes = new ArrayList<SuperTarget>();
			double distcovered = -1;
			if(a1.getNeighbors().contains(t))
			{
				//pnodes = base.getPath(dest);
				int src = a1.getTargetid();
				int des = t.getTargetid();

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];


			}
			else
			{




				int src = a1.getTargetid();
				int des = t.getTargetid();


				

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];
				//System.out.print("dist covered "+ distcovered+"\n");



				ArrayList<Integer>	tmppathnodes = apsp.getPath(src, des, apspmap, apspmapback);

				for(int k=0; k<tmppathnodes.size(); k++)
				{
					pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
				}

				

				//throw new Exception("Base to not neighbor for initial set of paths **********8");

			}
			
			
			if(dmin>distcovered)
			{

				dmin = distcovered;
				tid = t.getTargetid();
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				//tmppath.add(a1.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
				//	System.out.print(dest.getTargetid()+"\n");
				//tmppath.add(t.getTargetid());
				//System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				/*System.out.print("complete path : \n");
				for(int k=0; k<tmppath.size(); k++)
				{

					System.out.print(tmppath.get(k)+"->");
				}

				System.out.print("\n");
*/				
				minpath.clear();
				for(Integer p: tmppath)
				{
					minpath.add(p);
				}
				/*if(t.getTargetid()==0) //priority
				{
					break;
				}*/
			}

		}
		
		
		
		tmpa1path.clear();
		for(Integer p: minpath)
		{
			tmpa1path.add(p);
		}
		
		
		
		return new double[]{dmin, tid};
	}
	
	
	private static ArrayList<Integer>[] clusterWithWeka(int k, Instances newinstance, int totalcluster,
			ArrayList<Integer> targetstocluster, HashMap<Integer, TargetNode> targetmaps, EM dc,
			HashMap<Integer,Integer> apspmap, int[][] apspmat) throws Exception {
		
		
		
		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[totalcluster];
		
		for(int i=0; i<totalcluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}
		
		
	/*	FileInputStream fstream = new FileInputStream("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/realdata2.csv");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		Instances instances = new Instances(br);*/
		 
		 
		 
		 // Print header and instances.
		// System.out.println("\nDataset:\n");
		 //System.out.println(instances.toSummaryString());
		 
		/* SimpleKMeans model = new SimpleKMeans();
		 model.setNumClusters(10);
		 model.buildClusterer(instances);
		 model.setDistanceFunction(new weka.core.ManhattanDistance());
		 
		 System.out.println(model);*/
		
		
		//EM dc = new EM();
		
		 
		//SimpleKMeans dc = new SimpleKMeans();
		
		 
		 if(newinstance.get(0).value(0) != 0)
		 {
			 throw new Exception("0 is not the base");
		 }
		 
		 newinstance.remove(0); // remove base
		 
		 dc.setNumClusters(k-1); // 0 for base
		
		 dc.buildClusterer(newinstance);
		 System.out.println(dc);
		 
		 
		 
		 for(int i=0; i<newinstance.size(); i++)
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(newinstance.get(i))+1);
			 int clusterid = dc.clusterInstance(newinstance.get(i));
			 int tid = (int)newinstance.get(i).value(0);
			 if(tid!=0)
			 {
				 clusters[clusterid+1].add(tid);
			 }
		 
			 
		 }
		 
		 clusters[0].add(0); //base
		 
		// printClusters(clusters);
		 
		 int j = k;
		 for(Integer t: targetmaps.keySet())
		 {
			 if(!targetstocluster.contains(t))
			 {
				 clusters[j++].add(t.intValue());
			 }
		 }
		 

		// printClusters(clusters);
		 
		 // handle empty cluster
		 
		 
		 
		 
		
		return clusters;
	}

	
	

	private static void printClusters(ArrayList<Integer>[] cluster) {


		for(int i=0; i<cluster.length; i++)
		{
			System.out.print("cluster "+ i + ": ");
			for(Integer t: cluster[i])
			{
				System.out.print(t+ ", ");
			}
			System.out.println();
		}

	}
	
	
	private static void printClusters(HashMap<Integer, ArrayList<Integer>> clusters) {


		int i = 0;
		
		for(ArrayList<Integer> cluster: clusters.values())
		{
			System.out.print("cluster "+ i + ": ");
			for(Integer t: cluster)
			{
				System.out.print(t+ ", ");
			}
			i++;
			System.out.println();
		}

	}
	
	

	public static double[] wekaClusteringWithDO(int base, int dest, int ncluster, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, int iter, int slavelimit, int pathlimit) throws Exception
	{
		
		
		
		
		
		
		
		
		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		SecurityGameContraction.printSortedTargets(targetssorted);
		
		
		//Get the list of initial targets using GCR from Tsrt, Tcur = GreedyCoverR()
		
		HashMap<Integer, Integer> apspmap = new HashMap<>();
		HashMap<Integer, Integer> apspmapback = new HashMap<>();
		int[][] apspmat = new int[nTargets+1][nTargets+1];
		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		
		ArrayList<Integer> targetstocluster = buildGreedyCoverMultResWeka(targets, dmax, nTargets, 0, nRes,apspmap,
				apspmapback,apspmat, apsp );
		
		int currentPlace = targetstocluster.size()-1;
		
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		HashMap<Integer, TargetNode> tmptargetmaps = new HashMap<Integer, TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;

		


		long clusteringtime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];
		boolean canaddpath = true;
		
		int masteritr=0;
		
		int totalslaveiter = 0;
		
		
		
		
		/**
		 * make instances for clustering with weka
		 */
		
		
		
		
		

		 CSVLoader csvload = new CSVLoader();
		 csvload.setSource(new File("result/realdata"+iter+".csv"));
		 Instances data = csvload.getDataSet();
		 
		 
		 ArffSaver arf = new ArffSaver();
		 arf.setInstances(data);
		 
		 File f = new File("result/newdata"+iter+".arff");
		 
		 if(f.exists())
		 {
			 f.delete();
			 f.createNewFile();
		 }
		 
		 arf.setFile(f);
		 arf.writeBatch();
		
		
		
		
		
		
		
		FileInputStream fstream = new FileInputStream("result/newdata"+iter+".arff");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		// Read all the instances in the file (ARFF, CSV, XRFF, ...)
		 //DataSource source = new DataSource(br);
		 Instances instances = new Instances(br);
		 // Print header and instances.
		// System.out.println("\nDataset:\n");
		 System.out.println(instances.toSummaryString());
		 
		/* SimpleKMeans model = new SimpleKMeans();
		 model.setNumClusters(10);
		 model.buildClusterer(instances);
		 model.setDistanceFunction(new weka.core.ManhattanDistance());
		 System.out.println(model);*/
		 //MakeDensityBasedClusterer dc = new MakeDensityBasedClusterer();
		 EM dc = new EM();
		// instances.remove(0); // remove the base
		 dc.setNumClusters(ncluster-1);
		 dc.buildClusterer(instances);
		 System.out.println(dc);
		 
		/* for(int i=0; i<instances.size(); i++) //without base
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(instances.get(i)));
		  
		 }
		*/
		

		
		while(true)
		{
			ncluster = targetstocluster.size()/5;
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<targetstocluster.size(); i++)
			{
				System.out.print(targetstocluster.get(i)+",");
			}


			tmpgraph = SecurityGameContraction.getDuplicateGraph(targets);
			for(TargetNode t: tmpgraph)
			{
				tmptargetmaps.put(t.getTargetid(), t);
			}
			
			
			
			
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			//Build an abstract graph Gt through target clustering given Tcur, Gt
			
			
			HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
			HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
			
			
			if(!targetstocluster.contains(0))
			{
				System.out.println("No base mmmmmmm  ");

			}

			
			
			HashMap<Integer, SuperTarget> currentst = clusterTargetsWeka(targetstocluster, tmpgraph, 
					tmptargetmaps, dmax, ncluster, radius, dstravel, stpaths, dc, instances, apspmap, apspmat, apsp, apspmapback );
			
			targetsize= currentst.size();
			
			
			//printSuperTargets(currentst);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			clusteringtime += diff;
			
			
			//TODO save distance traveled for each cluster
			// remove the unncessary ones. 
			
			ArrayList<Integer> notin = new ArrayList<Integer>();
			
			for(SuperTarget st: currentst.values())
			{
				if(st.nodes.size()==1)
				{
					dstravel.put(st.stid, 0.0);
					
				}
				
			}
			
			
			int ind = 0;
			for(Integer t: dstravel.keySet())
			{
				if(!currentst.keySet().contains(t))
				{
					notin.add(t);
				}
				ind++;
			}
			
			for(Integer x: notin)
			{
				dstravel.remove(x);
				stpaths.remove(x);
			}
			
			
			//HashMap<Integer, Double> stvalue
			assignSTValues(currentst, tmptargetmaps);
			
			//System.out.println("olaa ");
			//printSuperTargets(currentst);
			
			
			//Generate initial set of paths using GreedyPathR, scur = GPR(Gt)
			
			

			

			System.out.println("current st size "+ currentst.size());
			
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[currentst.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, currentst, tmptargetmaps, nRes, dstravel);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: currentst.values())
			{

				map.put(st.stid, icount);
				System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//SecurityGameContraction.printPaths(pathseq);
			
			//System.out.println("hi");

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/
			ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();


			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;
				
				totalslaveiter++;

				

				if(pathseq.size()==0)
				{
					// might be because of no ap from ST 0. fix it
					System.out.println("pathseq 0, iter.ohhh"+ masteritr);
				}


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: currentst.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, currentst, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, currentst, tmptargetmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, currentst, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					
					System.out.println("attack target before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, currentst, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmptargetmaps, currentst, stpaths);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmptargetmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmptargetmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;
					
					int hiddenattackedclsuter = findAttackedCluster(currentst,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n  masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}


					
					

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();

					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmptargetmaps, 
							dmax, currentst.size(), 0, nRes, attackerstrategy, currentst, dstravel, pathlimit);
					
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, currentst, tmptargetmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);

					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
						System.out.println("newpathseq size before purify : "+newpathseq.size());
					    newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
						System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=slavelimit))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("newpathseq: ");
						//SecurityGameContraction.printPaths(newpathseq);

						System.out.println("Old path seq size "+ pathseq.size());

						int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						System.out.println("new paths added by slave *************");
						
						//SecurityGameContraction.printPaths(pathseq);

						/*pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						int newsize = pathseq.size();
						//System.out.println("haa ");


						if((oldsize==newsize) || (itr>=20))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ############### or iteration>20");
							printSuperTargets(currentst);
							SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				printSuperTargets(currentst);
				break;
			}
			
			double ulimit = tmptargetmaps.get(attackedtarget).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;
			int ADD_C = 5;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					targetstocluster.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=ADD_C)
					{
						break;
					}
				}
			}
			
			System.out.println("addcount : "+ addcount);

			currentPlace = targetstocluster.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<ADD_C || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += ADD_C-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					targetstocluster.add(targetssorted[k][0]);
				}
			}

			masteritr++;
			
		


		} // outer while loop

		System.out.println("Final target list size : "+ targetstocluster.size());

		for(int i=0; i<targetstocluster.size(); i++)
		{
			System.out.print(targetstocluster.get(i)+",");
		}

		
		double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, tmptargetmaps, probdistribution);


		//printSuperTargets(currentst);

		

		double[] res1 = {defpayoff, clusteringtime, solvingtime, currentPlace+1, attackeru, slavetime, revmaptime, totalslaveiter};
		return res1;
	}
	
	
	public static double[] wekaClusteringWithSO(int base, int dest, int ncluster, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, int iter, int nrow, int ncol, int slavelimit, int pathlimit) throws Exception
	{
		
		
		
		
		
		
		
		
		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		SecurityGameContraction.printSortedTargets(targetssorted);
		
		
		//Get the list of initial targets using GCR from Tsrt, Tcur = GreedyCoverR()
		
		HashMap<Integer, Integer> apspmap = new HashMap<>();
		HashMap<Integer, Integer> apspmapback = new HashMap<>();
		int[][] apspmat = new int[nTargets+1][nTargets+1];
		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		
		buildAPSPWeka(targets, dmax, nTargets, 0, nRes,apspmap,apspmapback,apspmat, apsp );
		
		
		ArrayList<Integer> targetstocluster = new ArrayList<Integer>();
		
		// include all the targets
		
		for(TargetNode t: targets)
		{
			targetstocluster.add(t.getTargetid());
		}
		
		
		int currentPlace = targetstocluster.size()-1;
		
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		HashMap<Integer, TargetNode> tmptargetmaps = new HashMap<Integer, TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;

		
		ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();

		long clusteringtime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];
		boolean canaddpath = true;
		
		int masteritr=0;
		
		int totalslaveiter = 0;
		
		
		
		
		/**
		 * make instances for clustering with weka
		 */
		
		
		
		
		

		 CSVLoader csvload = new CSVLoader();
		 csvload.setSource(new File("result/realdata"+iter+".csv"));
		 Instances data = csvload.getDataSet();
		 
		 
		 ArffSaver arf = new ArffSaver();
		 arf.setInstances(data);
		 
		 File f = new File("result/newdata"+iter+".arff");
		 
		 if(f.exists())
		 {
			 f.delete();
			 f.createNewFile();
		 }
		 
		 arf.setFile(f);
		 arf.writeBatch();
		
		
		
		
		
		
		
		FileInputStream fstream = new FileInputStream("result/newdata"+iter+".arff");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		// Read all the instances in the file (ARFF, CSV, XRFF, ...)
		 //DataSource source = new DataSource(br);
		 Instances instances = new Instances(br);
		 // Print header and instances.
		// System.out.println("\nDataset:\n");
		 System.out.println(instances.toSummaryString());
		 
		/* SimpleKMeans model = new SimpleKMeans();
		 model.setNumClusters(10);
		 model.buildClusterer(instances);
		 model.setDistanceFunction(new weka.core.ManhattanDistance());
		 System.out.println(model);*/
		 //MakeDensityBasedClusterer dc = new MakeDensityBasedClusterer();
		 EM dc = new EM();
		// instances.remove(0); // remove the base
		 dc.setNumClusters(ncluster-1);
		 dc.buildClusterer(instances);
		 System.out.println(dc);
		 
		/* for(int i=0; i<instances.size(); i++) //without base
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(instances.get(i)));
		  
		 }
		*/
		
		 HashMap<Integer, SuperTarget> currentst = new HashMap<Integer, SuperTarget>();
		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<targetstocluster.size(); i++)
			{
				System.out.print(targetstocluster.get(i)+",");
			}


			tmpgraph = SecurityGameContraction.getDuplicateGraph(targets);
			for(TargetNode t: tmpgraph)
			{
				tmptargetmaps.put(t.getTargetid(), t);
			}
			
			
			
			
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			//Build an abstract graph Gt through target clustering given Tcur, Gt
			
			
			HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
			HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
			
			
			if(!targetstocluster.contains(0))
			{
				System.out.println("No base mmmmmmm  ");

			}

			
			
			currentst = clusterTargetsWeka(targetstocluster, tmpgraph, 
					tmptargetmaps, dmax, ncluster, radius, dstravel, stpaths, dc, instances, apspmap, apspmat, apsp, apspmapback );
			
			targetsize= currentst.size();
			
			
			

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			clusteringtime += diff;
			
			
			//TODO save distance traveled for each cluster
			// remove the unncessary ones. 
			
			ArrayList<Integer> notin = new ArrayList<Integer>();
			
			for(SuperTarget st: currentst.values())
			{
				if(st.nodes.size()==1)
				{
					dstravel.put(st.stid, 0.0);
					
				}
				
			}
			
			
			int ind = 0;
			for(Integer t: dstravel.keySet())
			{
				if(!currentst.keySet().contains(t))
				{
					notin.add(t);
				}
				ind++;
			}
			
			for(Integer x: notin)
			{
				dstravel.remove(x);
				stpaths.remove(x);
			}
			
			
			//HashMap<Integer, Double> stvalue
			printSuperTargets(currentst, stpaths, dstravel);
			assignSTValues(currentst, tmptargetmaps);
			
			//System.out.println("olaa ");
			//printSuperTargets(currentst);
			
			
			//Generate initial set of paths using GreedyPathR, scur = GPR(Gt)
			
			

			

			System.out.println("current st size "+ currentst.size());
			
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[currentst.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, currentst, tmptargetmaps, nRes, dstravel);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: currentst.values())
			{

				map.put(st.stid, icount);
				System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




		//	SecurityGameContraction.printPaths(pathseq);
			
			
			if(pathseq.size()==0)
			{
			
				System.out.println("hi");
			}

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/

			currentattackedsupertargets.clear();

			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;
				
				totalslaveiter++;

				

				if(pathseq.size()==0)
				{
					// might be because of no ap from ST 0. fix it
					System.out.println("pathseq 0, iter.ohhh"+ masteritr);
				}


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: currentst.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, currentst, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, currentst, tmptargetmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, currentst, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					
					System.out.println("attack target before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, currentst, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmptargetmaps, currentst, stpaths);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmptargetmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmptargetmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;
					
					int hiddenattackedclsuter = findAttackedCluster(currentst,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n  masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}


					
					

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();

					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmptargetmaps, 
							dmax, currentst.size(), 0, nRes, attackerstrategy, currentst, dstravel, pathlimit);
					
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, currentst, tmptargetmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);

					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
						System.out.println("newpathseq size before purify : "+newpathseq.size());
					    newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
						System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=slavelimit))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("newpathseq: ");
						SecurityGameContraction.printPaths(newpathseq);

						System.out.println("Old path seq size "+ pathseq.size());

						int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						System.out.println("new paths added by slave *************");
						
						//SecurityGameContraction.printPaths(pathseq);

						/*pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						int newsize = pathseq.size();
						//System.out.println("haa ");


						if((oldsize==newsize) || (itr>=20))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ############### or iteration>20");
							printSuperTargets(currentst);
							SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				printSuperTargets(currentst);
				break;
			}
			
			double ulimit = tmptargetmaps.get(attackedtarget).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;
			int ADD_C = 5;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					targetstocluster.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=ADD_C)
					{
						break;
					}
				}
			}
			
			System.out.println("addcount : "+ addcount);

			currentPlace = targetstocluster.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<ADD_C || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += ADD_C-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					targetstocluster.add(targetssorted[k][0]);
				}
			}

			masteritr++;
			
		


		} // outer while loop

		System.out.println("Final target list size : "+ targetstocluster.size());

		for(int i=0; i<targetstocluster.size(); i++)
		{
			System.out.print(targetstocluster.get(i)+",");
		}

		
		double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, tmptargetmaps, probdistribution);


		//ButtonGrid grid = new ButtonGrid(targetmaps, currentst, "WekaWithSO");
		//grid.drawPayoffGrid(nrow, ncol);
		//grid.drawCluster(nrow, ncol);

		printSuperTargets(currentst);

		double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime, totalslaveiter};
		return res1;
	}
	
	
	public static double[] naiveClusteringWithSO(int base, int dest, int ncluster, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, int iter, int nrow, int ncol, int blocksize) throws Exception
	{
		
		
		
		
		
		
		
		
		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		SecurityGameContraction.printSortedTargets(targetssorted);
		
		
		//Get the list of initial targets using GCR from Tsrt, Tcur = GreedyCoverR()
		
		HashMap<Integer, Integer> apspmap = new HashMap<>();
		HashMap<Integer, Integer> apspmapback = new HashMap<>();
		int[][] apspmat = new int[nTargets+1][nTargets+1];
		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		
		buildAPSPWeka(targets, dmax, nTargets, 0, nRes,apspmap,apspmapback,apspmat, apsp );
		
		
		ArrayList<Integer> targetstocluster = new ArrayList<Integer>();
		
		// include all the targets
		
		for(TargetNode t: targets)
		{
			targetstocluster.add(t.getTargetid());
		}
		
		
		int currentPlace = targetstocluster.size()-1;
		
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		HashMap<Integer, TargetNode> tmptargetmaps = new HashMap<Integer, TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;

		
		ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();

		long clusteringtime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];
		boolean canaddpath = true;
		
		int masteritr=0;
		
		int totalslaveiter = 0;
		
		HashMap<Integer, SuperTarget> currentst = new HashMap<Integer, SuperTarget>();

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<targetstocluster.size(); i++)
			{
				System.out.print(targetstocluster.get(i)+",");
			}


			tmpgraph = SecurityGameContraction.getDuplicateGraph(targets);
			for(TargetNode t: tmpgraph)
			{
				tmptargetmaps.put(t.getTargetid(), t);
			}
			
			
			
			
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			//Build an abstract graph Gt through target clustering given Tcur, Gt
			
			
			HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
			HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
			
			
			if(!targetstocluster.contains(0))
			{
				System.out.println("No base mmmmmmm  ");

			}

			
			
			currentst = clusterTargetsNaive(targetstocluster, tmpgraph, 
					tmptargetmaps, dmax, ncluster, radius, dstravel, stpaths, apspmap, apspmat, apsp, apspmapback, nrow, blocksize );
			
			targetsize= currentst.size();
			
			
			//printSuperTargets(currentst);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			clusteringtime += diff;
			
			
			//TODO save distance traveled for each cluster
			// remove the unncessary ones. 
			
			ArrayList<Integer> notin = new ArrayList<Integer>();
			
			for(SuperTarget st: currentst.values())
			{
				if(st.nodes.size()==1)
				{
					dstravel.put(st.stid, 0.0);
					
				}
				
			}
			
			
			int ind = 0;
			for(Integer t: dstravel.keySet())
			{
				if(!currentst.keySet().contains(t))
				{
					notin.add(t);
				}
				ind++;
			}
			
			for(Integer x: notin)
			{
				dstravel.remove(x);
				stpaths.remove(x);
			}
			
			
			//HashMap<Integer, Double> stvalue
			assignSTValues(currentst, tmptargetmaps);
			
			//System.out.println("olaa ");
			//printSuperTargets(currentst);
			
			
			//Generate initial set of paths using GreedyPathR, scur = GPR(Gt)
			
			

			

			System.out.println("current st size "+ currentst.size());
			
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[currentst.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, currentst, tmptargetmaps, nRes, dstravel);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: currentst.values())
			{

				map.put(st.stid, icount);
				System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




		//	SecurityGameContraction.printPaths(pathseq);
			
			
			if(pathseq.size()==0)
			{
			
				System.out.println("hi");
			}

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/

			currentattackedsupertargets.clear();

			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;
				
				totalslaveiter++;

				

				if(pathseq.size()==0)
				{
					// might be because of no ap from ST 0. fix it
					System.out.println("pathseq 0, iter.ohhh"+ masteritr);
				}


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: currentst.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, currentst, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, currentst, tmptargetmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, currentst, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					
					System.out.println("attack target before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, currentst, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmptargetmaps, currentst, stpaths);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmptargetmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmptargetmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;
					
					int hiddenattackedclsuter = findAttackedCluster(currentst,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n  masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}


					
					

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();

					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSuperGreedyCoverMultRes2(tmptargetmaps, 
							dmax, currentst.size(), 0, nRes, attackerstrategy, currentst, dstravel);
					
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, currentst, tmptargetmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);

					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
						System.out.println("newpathseq size before purify : "+newpathseq.size());
					    newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
						System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=10))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("newpathseq: ");
						SecurityGameContraction.printPaths(newpathseq);

						System.out.println("Old path seq size "+ pathseq.size());

						int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						System.out.println("new paths added by slave *************");
						
						//SecurityGameContraction.printPaths(pathseq);

						/*pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						int newsize = pathseq.size();
						//System.out.println("haa ");


						if((oldsize==newsize) || (itr>=20))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ############### or iteration>20");
							printSuperTargets(currentst);
							SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				printSuperTargets(currentst);
				break;
			}
			
			double ulimit = tmptargetmaps.get(attackedtarget).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;
			int ADD_C = 5;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					targetstocluster.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=ADD_C)
					{
						break;
					}
				}
			}
			
			System.out.println("addcount : "+ addcount);

			currentPlace = targetstocluster.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<ADD_C || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += ADD_C-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					targetstocluster.add(targetssorted[k][0]);
				}
			}

			masteritr++;
			
		


		} // outer while loop

		System.out.println("Final target list size : "+ targetstocluster.size());

		for(int i=0; i<targetstocluster.size(); i++)
		{
			System.out.print(targetstocluster.get(i)+",");
		}

		
		double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, tmptargetmaps, probdistribution);


		ButtonGrid grid = new ButtonGrid(targetmaps, currentst, "NaiveWithSO");
		//grid.drawPayoffGrid(nrow, ncol);
		grid.drawCluster(nrow, ncol);

		

		double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime, totalslaveiter};
		return res1;
	}


	public static int[][] makeSuperOrigPMatWOMap(int[][] p,
			ArrayList<ArrayList<Integer>> pathseq, 
			List<ArrayList<Integer>> jset, int nTargets,
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback, 
			HashMap<Integer,TargetNode> targetmaps, HashMap<Integer,SuperTarget> currentst, HashMap<Integer,ArrayList<Integer>> stpaths) throws Exception {



		int[][] origpmat = new int[nTargets][jset.size()];

		int jindex = 0;
		for(ArrayList<Integer> j: jset)
		{

			ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
			for(Integer jpath: j)
			{
				paths.add(pathseq.get(jpath));
			}


			// 0 11 18 2 5 9 13 19 6 7 12 14 15 16 20 21 22 24 

			for(ArrayList<Integer> path: paths)
			{
				//printNodesWithNeighborsAndPath(dominatednodes, targets);
				/*System.out.println("Considering path :");
				printGreedyPath(path);
				System.out.println();*/

				for(int targetindex=0; targetindex<path.size()-1; targetindex++)
				{

					// get the supertarget
					int st = path.get(targetindex);
					int nextst = path.get(targetindex+1);
					
					// for all the nodes in the supertarget assign 1
					
					for(TargetNode t: currentst.get(st).nodes.values())
					{
						
						origpmat[t.getTargetid()][jindex] = 1;
						
							
					}
					
					for(TargetNode t: currentst.get(nextst).nodes.values())
					{
						origpmat[t.getTargetid()][jindex] = 1;
					}
					
					
					// get the shortest path between st and nextst
					// then cover those path nodes between two ap
					
					ArrayList<TargetNode> inbetweennodes = inBetweenNodes(st,nextst, currentst);
					
					for(TargetNode n: inbetweennodes)
					{
						origpmat[n.getTargetid()][jindex] = 1;
					}
					
					
				}

			}
			jindex++;
		}


		return origpmat;
	}


	private static ArrayList<TargetNode> inBetweenNodes(int st, int nextst, HashMap<Integer,SuperTarget> currentst) {
		
		
		ArrayList<TargetNode> nodes = new ArrayList<TargetNode>();
		
		
		SuperTarget t1 = currentst.get(st);
		SuperTarget t2 = currentst.get(nextst);
		
		int t1ap = -1;
		int t2ap= -1;
		double mindist = Double.MAX_VALUE;
		
		for(TargetNode n1: t1.ap.values())
		{
			for(TargetNode n2: t2.ap.values())
			{
				
				if(n1.getNeighbors().contains(n2))
				{
					double dist = n1.getDistance(n2);
					
					if(dist<mindist)
					{
						mindist = dist;
						nodes = n1.getPath(n2);
					}
				}
			}
		}
		
		
		return nodes;
		
		
		
	}
	
	
public static ArrayList<TargetNode> inBetweenNodes(SuperTarget t1, SuperTarget t2) {
		
		
		ArrayList<TargetNode> nodes = new ArrayList<TargetNode>();
		
		
		
		
		int t1ap = -1;
		int t2ap= -1;
		double mindist = Double.MAX_VALUE;
		
		for(TargetNode n1: t1.ap.values())
		{
			for(TargetNode n2: t2.ap.values())
			{
				
				if(n1.getNeighbors().contains(n2))
				{
					double dist = n1.getDistance(n2);
					
					if(dist<mindist)
					{
						mindist = dist;
						nodes = n1.getPath(n2);
					}
				}
			}
		}
		
		
		return nodes;
		
		
		
	}


	private static boolean isDone(ArrayList<int[]> done, TargetNode n, TargetNode n2) {
		
		
		for(int[] x: done)
		{
			if((x[0] == n.getTargetid() && x[1] == n2.getTargetid()) ||(x[1] == n.getTargetid() && x[0] == n2.getTargetid()))
				return true;
		}
		
		
		return false;
	}


	

	

	private static ArrayList<Integer> makeAttackCluster(int base, int attackedtarget,
			HashMap<Integer, TargetNode> targetmaps) {
		// TODO Auto-generated method stub
		
		
		ArrayList<Integer> clus = new ArrayList<Integer>();
		
		
		for(TargetNode t: targetmaps.get(attackedtarget).getNeighbors())
		{
			if(t.getTargetid() != base)
			{
				t.isinattackcluster = true;
				clus.add(t.getTargetid());
			}
		}
		
		
		return clus;
		
		
	}
	
	
	
public static int[][] constructGameData(ArrayList<TargetNode> u) {
		
		
		int[][] gd = new int[u.size()][u.size()];
		
		int index =0;
		
		for(int i=0; i<u.size(); i++)
		{
			//for(int j=0; j<u.size(); j++)
			{
				gd[index][0] = 0;
				gd[index][1] = (int)Math.floor(-u.get(i).attackerreward);
				gd[index][2] = -gd[index][1];
				gd[index][3] = 0;
				index++;
			}
		}
		
		
		
		return gd;
	}
	

public static void generateGraphWithOutContraction(ArrayList<TargetNode> tmpgraph, int currentPlace, ArrayList<TargetNode> domindatednodes, 
		int[][] targetssorted, double dmax, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> tmpgraphmaps)
{

	
	ArrayList<TargetNode> graph = SecurityGameContraction.getDuplicateGraph(targets);
	
	
	for(TargetNode t: graph)
	{
		tmpgraphmaps.put(t.getTargetid(), t);
	}
	
	
	ArrayList<TargetNode> dom = new ArrayList<TargetNode>();
	if(currentPlace<targetssorted.length-1)
		dom = SecurityGameContraction.selectDominatedNodes(targetssorted, currentPlace+1, graph);
	else
	{
		dom.clear();
	}

//	System.out.print("\nDom targets : ");
	for(TargetNode s: dom)
	{
		//System.out.print(s.getTargetid()+" ");
		domindatednodes.add(s);
	}
	//System.out.println();

	//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


	//System.out.println("Contracting graph...");

	//SecurityGameContraction.instantContractionWithAPSP(domindatednodes, graph, dmax);
	
	//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

	//SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, graph);
	SecurityGameContraction.removeDominatedTargets(domindatednodes, graph);
	
	for(TargetNode t: graph)
	{
		tmpgraph.add(t);
	}
	
	
	
	//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
	

	System.out.println("tmpgraph size "+ tmpgraph.size());
	System.out.println("dom size "+ domindatednodes.size());
	//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

	
	}


public static void generateGraph(ArrayList<TargetNode> tmpgraph, int currentPlace, ArrayList<TargetNode> domindatednodes, 
		int[][] targetssorted, double dmax, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> tmpgraphmaps)
{

	
	ArrayList<TargetNode> graph = SecurityGameContraction.getDuplicateGraph(targets);
	
	
	for(TargetNode t: graph)
	{
		tmpgraphmaps.put(t.getTargetid(), t);
	}
	
	
	ArrayList<TargetNode> dom = new ArrayList<TargetNode>();
	if(currentPlace<targetssorted.length-1)
		dom = SecurityGameContraction.selectDominatedNodes(targetssorted, currentPlace+1, graph);
	else
	{
		dom.clear();
	}

	System.out.print("\nDom targets : ");
	for(TargetNode s: dom)
	{
		//System.out.print(s.getTargetid()+" ");
		domindatednodes.add(s);
	}
	System.out.println();

	//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


	System.out.println("Contracting graph...");

	SecurityGameContraction.instantContractionWithAPSP(domindatednodes, graph, dmax);
	
	//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

	SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, graph);
	SecurityGameContraction.removeDominatedTargets(domindatednodes, graph);
	
	for(TargetNode t: graph)
	{
		tmpgraph.add(t);
	}
	
	
	
	//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
	

	System.out.println("tmpgraph size "+ tmpgraph.size());
	System.out.println("dom size "+ domindatednodes.size());
	//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

	
	}



public static ArrayList<ArrayList<Integer>>  generatePaths(ArrayList<TargetNode> tmpgraph, ArrayList<TargetNode> targets,
		double dmax, ArrayList<Integer> currenttargets, int currentPlace, ArrayList<ArrayList<Integer>> pathsequence,
		int[][] gamedata, HashMap<Integer, Integer> mapback, HashMap<Integer, Integer> map, int nRes) throws Exception
{

	
	
	//apply greedy approach
	//TODO generate paths where there will be at least one target
	//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
	//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
	
	
	
	ArrayList<ArrayList<Integer>> pathseq =  SecurityGameContraction.generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
	//map = new HashMap<Integer, Integer>();
	//mapback = new HashMap<Integer, Integer>();
	int icount =0;
	for(int i=0; i<tmpgraph.size(); i++)
	{

		map.put(tmpgraph.get(i).getTargetid(), icount);
		//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
		mapback.put(icount, tmpgraph.get(i).getTargetid());
		icount++;

	}
	//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
	//SecurityGameContraction.printPaths(pathseq);
	System.out.println("Total path with duplicates "+pathseq.size());
	pathseq =SecurityGameContraction.removeDuplicatePathSimple(pathseq);
	System.out.println("Total path without duplicates "+pathseq.size()+"\n");
	SecurityGameContraction.printPaths(pathseq);
	
	

	/**
	 * keep only nRes*3 paths from the end
	 */

	//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
	//System.out.println("Initial number of paths "+ pathseq.size());
	//printPaths(pathseq);


	/*if(pathseq.size()==0)
	{
		System.out.println("pathseq 0, iter..nqqqqq"+ iter);
	}*/
	
	return pathseq;


	
	}




public static double pathAPForST( HashMap<Integer,TargetNode> nodes, double dmax, int nTargets, int base, int dest, int nRes, ArrayList<Integer> tmpa1a2spath) {



	int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	for(TargetNode t: nodes.values())
	{
		targets.add(t);
	}
	

	/**
	 * make mapping
	 */

	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	int icount =1;
	for(int i=0; i<targets.size(); i++)
	{
		map.put(targets.get(i).getTargetid(), icount);
		//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
		mapback.put(icount, targets.get(i).getTargetid());
		icount++;
	}
	//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

	makeAdjacencyMatrixST(adjacencymatrix, targets, nTargets, map, mapback, nodes);


	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
	
	
	SecurityGameContraction.purifyAPSPMatrixZero(apspmat, targets, nTargets, map, mapback);
	



	ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	
	
	int s = map.get(base);
	int d = map.get(dest);
	dmax = apspmat[s][d];
	//disallocation[0] = (int)dmax;

	tcur = SecurityGameContraction.greedyOnePathWithSrcDest(base, dest, targets, dmax, targetssorted, apspmat, map,mapback, nRes);
	
	/*for(Integer n: tcur)
	{
		tmpa1a2spath.add(n);
	}*/

	return tcur.size();
}





	

private static double[] DOWithPACMANClus(int[][] gamedata,
		int nTargets, int nRes, double[][] density, double
		dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int RADIUS, HashMap<Integer, Integer> clusterhistogram, int slavelimit, int pathlimit) throws Exception {


	

	/*targets.clear();
	SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
	
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

	//printtargets(targets);

	/**
	 * 1. sort the targets
	 */
	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	SecurityGameContraction.printSortedTargets(targetssorted);
	
	
	HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>();
	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	

	HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
	HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
	
	
	
	int apspmat[][] =  buildAPSP(targets, nTargets, apspmap, apspmapback, apsp);
	
	

	ArrayList<Integer> currenttargets = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
	//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
	/*currenttargets.add(targetssorted[0][0]);
	currenttargets.add(targetssorted[1][0]);*/



	int currentPlace = currenttargets.size()-1;


	ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

	ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
	int attackedtarget=-1;
	ArrayList<Integer> attackhistory = new ArrayList<Integer>();
	int[][] p;
	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	HashSet jSet=new HashSet();
	ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
	List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
	double[] probdistribution= new double[jset.size()];
	double attackeru= -999;
	double attackerv = -999;

	


	long contractiontime=0;
	long solvingtime=0;
	long revmaptime=0;
	int targetsize=0;
	double slavetime = 0;
	int [][] origpmat = new int[nTargets][];
	int totalslaveiter = 0;
	long clusteringtime = 0;






	boolean canaddpath = true;
	boolean clusteringactivated = false;
	int masteritr=0;

	
	
	HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();
	
	
	
	// cluster history which will be maintained always
	HashMap<Integer, ArrayList<Integer>> attackclustershisotry = new HashMap<Integer, ArrayList<Integer>>();
	
	// list of targets which has been attacked
	ArrayList<Integer> clusteredtargets = new ArrayList<Integer>();
	
	// new targets added for attacker strategy
	// start with the current targets
	// then only add new targets selected by attacker oracle
	ArrayList<Integer> newtargetstocluster = new ArrayList<Integer>();
	
	// first target is the center of a cluster
	HashMap<Integer,Integer> clustercenters = new HashMap<Integer,Integer>();
	
	
	//cluster access points 
	HashMap<Integer, int[]> clusterap = new HashMap<Integer, int[]>();
	
	
	// current attacked Super targets in defender oracle
	ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();
	
	
	for(Integer t: currenttargets)
	{
		newtargetstocluster.add(t);
	}
	
	
	//ArrayList<Integer> notaddedst = new ArrayList<Integer>(); // ST which were not added by slave
	
	
	
	
	while(true)
	{
		
		
		
		
		System.out.println("\n clusteringactivated "+clusteringactivated +" Outer loop...Master");

		pathseq = new ArrayList<ArrayList<Integer>>();

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current place : "+ currentPlace);

		System.out.print("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}
		
		p = new int[targets.size()][]; // p matrix



		Date start = new Date();
		long l1 = start.getTime();

		
		tmpgraph.clear();
		domindatednodes.clear();
		HashMap<Integer, TargetNode> tmpgraphmaps = new HashMap<Integer, TargetNode>();
		
		generateGraph(tmpgraph, currentPlace, domindatednodes, targetssorted, dmax, targets, tmpgraphmaps);
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps);
		
		
		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		
		map.clear();
		mapback.clear();
		
		clusteringactivated = true;// clsutering activated from the very beiginning

		
		
		
		if(clusteringactivated)
		{
			System.out.println("need clustering "+clusteringactivated+", graph size "+ tmpgraph.size());
			//  build a graph with super targets
			// build an attack cluster except the dominated targets
			//sts.clear();
			dstravel.clear();
			stpaths.clear();
			
			
			
			start = new Date();
			l1 = start.getTime();

			


			
			sts = constructSuperTargets(tmpgraphmaps, attackhistory, apsp, apspmat, apspmap,apspmapback, dstravel, stpaths, (int) dmax, attackclustershisotry,
					clusteredtargets, clustercenters, newtargetstocluster, domindatednodes, RADIUS, tmpgraph, clusterap);
			
			

			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			clusteringtime += diff;
			

			printSuperTargets(sts, stpaths, dstravel);
			preparePaths(dstravel, stpaths, sts);
			assignSTValues(sts, tmpgraphmaps);
			
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel);
			
			
			if(pathseq.size()==0)
			{
				System.out.println("no path seq");
			}
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths before removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths after removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: sts.values())
			{

				map.put(st.stid, icount);
				//System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			
			
			
			
			
			
		}

		

		int itr=0;
		
		
		newtargetstocluster.clear();
		
		currentattackedsupertargets.clear();
		
		//notaddedst.clear();
		
		
		while(true)
		{
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", Entered inner loop...slave");
			
			itr++;
			totalslaveiter++;

			

			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter.ohhh"+ iter);
			}*/


			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			
			
			// end if for clustering activated
			if(clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: sts.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, sts, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, sts, tmpgraphmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, sts, map, mapback);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target before rev map "+ attackedtarget);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack ST before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, sts, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmpgraphmaps, sts, stpaths);
					
					
					
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmpgraphmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmpgraphmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target after rev map"+ attackedtarget);
					
					
					/*if(!newtargetstocluster.contains(attackedtarget))
					{
						newtargetstocluster.add(attackedtarget);
					}*/
					
					
					
					// find the attacked cluster
					// add it to the attackedsupertarget
					
					int hiddenattackedclsuter = findAttackedCluster(sts,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}
					
					
					
					
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +", inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();
					
					//printSuperTargets(sts);
					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmpgraphmaps, 
							dmax, sts.size(), 0, nRes, attackerstrategy, sts, dstravel, pathlimit);
					
					
					
					//printSuperTargets(sts);
					
					
					//////////// ADD Attacked clusters in the path////////
					SecurityGameContraction.addSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);
					

					//////////////////////////////////////





					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;


					System.out.println("newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					
					/**
					 * check if the hidden attacked cluster or attacked clusters were added to newpathseq
					 */
					
					
					//notaddedst = getNotAddedST(newpathseq, currentattackedsupertargets);
					


					if( ((newpathseq.size()==0) || (itr>=pathlimit)))
					{
						canaddpath = false;
						System.out.println("not ST left to add to path....Slave can't add any new path ###############");
						break;
					}
					
					
					System.out.println("New whole path seq ");


					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Old path seq size "+ pathseq.size());

					
					
					
					
					
					//int oldsize = pathseq.size();

					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}
					
					//SecurityGameContraction.removeDuplicatePathSimple(pathseq);
					
					//int newsize = pathseq.size();
					
					/*if(oldsize == newsize)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" New path seq size "+ pathseq.size());


					SecurityGameContraction.printPaths(pathseq);

				} // end if else
			}
			System.out.println("iter"+ itr);
			
		} // inner while loop 
		
		newtargetstocluster.clear();// clear it to add new targets for clustering




		// add all targets all targets with utility >= U(a')


		if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
			break;
		}

		


		double ulimit = SecurityGameContraction.getTargetNode(attackedtarget, targets).attackerreward;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


		int addcount=0;

		for(int k=currentPlace+1; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				addcount++;
				
				currenttargets.add(targetssorted[k][0]);
				newtargetstocluster.add(targetssorted[k][0]);
				System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
				if(addcount>=5)
				{
					break;
				}
			}
		}

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " addcount : "+ addcount);

		currentPlace = currenttargets.size()-1;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " currentplace  : "+ currentPlace);

		if(addcount<5 || addcount==0)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding more ");

			int prevcur = currentPlace;
			currentPlace += 5-addcount;

			//System.out.println("currentplace  : "+ currentPlace);
			if(currentPlace>targetssorted.length-1)
			{
				currentPlace = targetssorted.length-1;
			}
			//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{

				System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding target  "+ targetssorted[k][0]);
				currenttargets.add(targetssorted[k][0]);
				newtargetstocluster.add(targetssorted[k][0]);
			}
		}

		masteritr++;


	} // outer while loop

	System.out.println("Final target list size : "+ currenttargets.size());

	for(int i=0; i<currenttargets.size(); i++)
	{
		System.out.print(currenttargets.get(i)+",");
	}

	//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
	double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);

	System.out.println("def exp : "+ defpayoff);
	
	
	/*printST(sts, nTargets, iter);
	
	printClusteredNodes(sts,nTargets,iter);
	
	printClusterDists(sts,nTargets,iter);
	*/
	
	for(SuperTarget st: sts.values())
	{
		int value = st.nodes.size();
		if(clusterhistogram.containsKey(st.stid))
		{
			value += clusterhistogram.get(st.stid);
			
		}
		clusterhistogram.put(st.stid, value);
	}
	
	
	
//	verifySolution(jset,pathseq,probdistribution, nTargets, dmax, sts, dstravel, tmpgraph);

	
	
	//ButtonGrid grid = new ButtonGrid(targetmaps, sts, "DOWithPACMAN");
	//grid.drawPayoffGrid(nrow, ncol);
	//grid.drawCluster(nrow, ncol);


	//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

	double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime, totalslaveiter, clusteringtime};
	return res;
}


private static double[] DOWithSplitPACMANClus(int[][] gamedata,
		int nTargets, int nRes, double[][] density, double
		dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int RADIUS, HashMap<Integer, Integer> clusterhistogram, int slavelimit, int pathlimit) throws Exception {


	

	/*targets.clear();
	SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
	
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

	//printtargets(targets);

	/**
	 * 1. sort the targets
	 */
	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	SecurityGameContraction.printSortedTargets(targetssorted);
	
	
	HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>();
	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	

	HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
	HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
	
	
	
	int apspmat[][] =  buildAPSP(targets, nTargets, apspmap, apspmapback, apsp);
	
	

	ArrayList<Integer> currenttargets = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
	//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
	/*currenttargets.add(targetssorted[0][0]);
	currenttargets.add(targetssorted[1][0]);*/



	int currentPlace = currenttargets.size()-1;


	ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

	ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
	int attackedtarget=-1;
	ArrayList<Integer> attackhistory = new ArrayList<Integer>();
	int[][] p;
	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	HashSet jSet=new HashSet();
	ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
	List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
	double[] probdistribution= new double[jset.size()];
	double attackeru= -999;
	double attackerv = -999;

	


	long contractiontime=0;
	long solvingtime=0;
	long revmaptime=0;
	int targetsize=0;
	double slavetime = 0;
	int [][] origpmat = new int[nTargets][];
	int totalslaveiter = 0;
	long clusteringtime = 0;






	boolean canaddpath = true;
	boolean clusteringactivated = false;
	int masteritr=0;

	
	
	HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();
	
	
	
	// cluster history which will be maintained always
	HashMap<Integer, ArrayList<Integer>> attackclustershisotry = new HashMap<Integer, ArrayList<Integer>>();
	
	// list of targets which has been attacked
	ArrayList<Integer> clusteredtargets = new ArrayList<Integer>();
	
	// new targets added for attacker strategy
	// start with the current targets
	// then only add new targets selected by attacker oracle
	ArrayList<Integer> newtargetstocluster = new ArrayList<Integer>();
	
	// first target is the center of a cluster
	HashMap<Integer,Integer> clustercenters = new HashMap<Integer,Integer>();
	
	
	//cluster access points 
	HashMap<Integer, int[]> clusterap = new HashMap<Integer, int[]>();
	
	
	// current attacked Super targets in defender oracle
	ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();
	
	
	for(Integer t: currenttargets)
	{
		newtargetstocluster.add(t);
	}
	
	
	ArrayList<Integer> notaddedst = new ArrayList<Integer>(); // ST which were not added by slave
	
	
	
	
	while(true)
	{
		
		
		
		
		System.out.println("\n clusteringactivated "+clusteringactivated +" Outer loop...Master");

		pathseq = new ArrayList<ArrayList<Integer>>();

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current place : "+ currentPlace);

		System.out.print("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}
		
		p = new int[targets.size()][]; // p matrix



		Date start = new Date();
		long l1 = start.getTime();

		
		tmpgraph.clear();
		domindatednodes.clear();
		HashMap<Integer, TargetNode> tmpgraphmaps = new HashMap<Integer, TargetNode>();
		
		generateGraph(tmpgraph, currentPlace, domindatednodes, targetssorted, dmax, targets, tmpgraphmaps);
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps);
		
		
		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		
		map.clear();
		mapback.clear();
		
		clusteringactivated = true;// clsutering activated from the very beiginning

		
		
		
		if(clusteringactivated)
		{
			System.out.println("need clustering "+clusteringactivated+", graph size "+ tmpgraph.size());
			//  build a graph with super targets
			// build an attack cluster except the dominated targets
			//sts.clear();
			dstravel.clear();
			stpaths.clear();
			
			
			
			start = new Date();
			l1 = start.getTime();

			


			
			sts = constructSuperTargetsV2(tmpgraphmaps, attackhistory, apsp, apspmat, apspmap,apspmapback, dstravel, stpaths, (int) dmax, attackclustershisotry,
					clusteredtargets, clustercenters, newtargetstocluster, domindatednodes, RADIUS, tmpgraph, clusterap, notaddedst);
			
			

			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			clusteringtime += diff;
			

			//printSuperTargets(sts, stpaths, dstravel);
			preparePaths(dstravel, stpaths, sts);
			assignSTValues(sts, tmpgraphmaps);
			
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel);
			
			
			if(pathseq.size()==0)
			{
				System.out.println("no path seq");
			}
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths before removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths after removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: sts.values())
			{

				map.put(st.stid, icount);
				//System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			
			
			
			
			
			
		}

		

		int itr=0;
		
		
		newtargetstocluster.clear();
		
		currentattackedsupertargets.clear();
		
		notaddedst.clear();
		
		
		while(true)
		{
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", Entered inner loop...slave");
			
			itr++;
			totalslaveiter++;

			

			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter.ohhh"+ iter);
			}*/


			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			
			
			// end if for clustering activated
			if(clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: sts.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, sts, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, sts, tmpgraphmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, sts, map, mapback);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target before rev map "+ attackedtarget);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack ST before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, sts, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmpgraphmaps, sts, stpaths);
					
					
					
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmpgraphmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmpgraphmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target after rev map"+ attackedtarget);
					
					
					/*if(!newtargetstocluster.contains(attackedtarget))
					{
						newtargetstocluster.add(attackedtarget);
					}*/
					
					
					
					// find the attacked cluster
					// add it to the attackedsupertarget
					
					int hiddenattackedclsuter = findAttackedCluster(sts,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}
					
					
					
					
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +", inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();
					
					//printSuperTargets(sts);
					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmpgraphmaps, 
							dmax, sts.size(), 0, nRes, attackerstrategy, sts, dstravel, pathlimit);
					
					
					
					printSuperTargets(sts);
					
					
					//////////// ADD Attacked clusters in the path////////
					SecurityGameContraction.addSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);
					

					//////////////////////////////////////





					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;


					System.out.println("newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					
					/**
					 * check if the hidden attacked cluster or attacked clusters were added to newpathseq
					 */
					
					
					
					

					

					/*if((notaddedst.size()==0) && ((newpathseq.size()==0) || (itr>=20)))
					{
						canaddpath = false;
						System.out.println("not ST left to add to path....Slave can't add any new path ###############");
						break;
					}*/
					
					
					if( ((newpathseq.size()==0) || (itr>=slavelimit)))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");


					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Old path seq size "+ pathseq.size());

					
					
					
					
					
					//int oldsize = pathseq.size();

					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}
					
					notaddedst = getNotAddedST(pathseq, currentattackedsupertargets);
					if(notaddedst.size()>0)
					{
						break;
					}
					
					//SecurityGameContraction.removeDuplicatePathSimple(pathseq);
					
					//int newsize = pathseq.size();
					
					/*if(oldsize == newsize)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" New path seq size "+ pathseq.size());


					SecurityGameContraction.printPaths(pathseq);

				} // end if else
			}
			System.out.println("iter"+ itr);
			
		} // inner while loop 
		
		newtargetstocluster.clear();// clear it to add new targets for clustering




		// add all targets all targets with utility >= U(a')


		if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
			break;
		}

		


		double ulimit = SecurityGameContraction.getTargetNode(attackedtarget, targets).attackerreward;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


		int addcount=0;

		for(int k=currentPlace+1; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				addcount++;
				
				currenttargets.add(targetssorted[k][0]);
				newtargetstocluster.add(targetssorted[k][0]);
				System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
				if(addcount>=5)
				{
					break;
				}
			}
		}

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " addcount : "+ addcount);

		currentPlace = currenttargets.size()-1;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " currentplace  : "+ currentPlace);

		if(addcount<5 || addcount==0)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding more ");

			int prevcur = currentPlace;
			currentPlace += 5-addcount;

			//System.out.println("currentplace  : "+ currentPlace);
			if(currentPlace>targetssorted.length-1)
			{
				currentPlace = targetssorted.length-1;
			}
			//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{

				System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding target  "+ targetssorted[k][0]);
				currenttargets.add(targetssorted[k][0]);
				newtargetstocluster.add(targetssorted[k][0]);
			}
		}

		masteritr++;


	} // outer while loop

	System.out.println("Final target list size : "+ currenttargets.size());

	for(int i=0; i<currenttargets.size(); i++)
	{
		System.out.print(currenttargets.get(i)+",");
	}

	//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
	double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);

	System.out.println("def exp : "+ defpayoff);
	
	
	/*printST(sts, nTargets, iter);
	
	printClusteredNodes(sts,nTargets,iter);
	
	printClusterDists(sts,nTargets,iter);
	*/
	
	for(SuperTarget st: sts.values())
	{
		int value = st.nodes.size();
		if(clusterhistogram.containsKey(st.stid))
		{
			value += clusterhistogram.get(st.stid);
			
		}
		clusterhistogram.put(st.stid, value);
	}
	
	
	
//	verifySolution(jset,pathseq,probdistribution, nTargets, dmax, sts, dstravel, tmpgraph);

	
	
	//ButtonGrid grid = new ButtonGrid(targetmaps, sts, "DOWithPACMAN");
	//grid.drawPayoffGrid(nrow, ncol);
	//grid.drawCluster(nrow, ncol);


	//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

	double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime, totalslaveiter, clusteringtime};
	return res;
}




private static ArrayList<Integer> getNotAddedST(ArrayList<ArrayList<Integer>> newpathseq,
		ArrayList<Integer> currentattackedsupertargets) {
	
	
	ArrayList<Integer> notadded = new ArrayList<Integer>();
	
	
	for(Integer st: currentattackedsupertargets)
	{
		boolean isthere = false;
		for(ArrayList<Integer> path: newpathseq)
		{
			if(path.contains(st))
			{
				isthere = true;
				break;
			}
		}
		if(!isthere)
			notadded.add(st);
	}
	
	
	
	return notadded;
}


public static void verifySolution(List<ArrayList<Integer>> jset, ArrayList<ArrayList<Integer>> pathseq,
		double[] probdistribution, int nTargets, double dmax, HashMap<Integer,SuperTarget> sts, HashMap<Integer,Double> dstravel, ArrayList<TargetNode> tmpgraph)
{
	
	
	


	ArrayList<Integer> donepaths = new ArrayList<Integer>();



	for(int probindex=0; probindex<probdistribution.length; probindex++)
	{
		if(probdistribution[probindex]>0)
		{
			for(int pathindex: jset.get(probindex))
			{


				if(!donepaths.contains(pathindex))
				{	

					donepaths.add(pathindex);

					System.out.print("\n\nPath "+ pathindex + ", prob: "+ probdistribution[probindex] + " ");




					ArrayList<Integer> path = pathseq.get(pathindex);

					double dist = 0;
					
					

					for(int index=0; index<path.size()-1; index++)
					{

						int pathnode1 = path.get(index);
						int pathnode2 = path.get(index+1);



						SuperTarget p1 =  sts.get(pathnode1);
						SuperTarget p2 =  sts.get(pathnode2);
						
						
						/*for(SuperTarget n: p1.distances.keySet())
						{
							System.out.println(n.stid+ " -> "+ p1.distances.get(n));
						}*/

						/*
						if(!p1.distances.keySet().contains(p2))
						{
							dist += getDist(p1,p2, tmpgraph);
						}
						else
						{
							dist +=  p1.distances.get(p2);
						}*/
						
						
						dist +=  p1.distances.get(p2);
						
						
						if(p2.nodes.size()>1)
						{
							dist += dstravel.get(p2.stid);
						}
						
						System.out.print(pathnode1+ "->");

					}


					System.out.println(" dist "+dist);
				}

			}
		}
	}

		
		//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
		//pw.append(expno+","+nTargets+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," +slavetime+","+ totaltime+"\n");
		//pw.close();

	
	
}


private static double getDist(SuperTarget p1, SuperTarget p2, ArrayList<TargetNode> tmpgraph) {
	
	
	Double min = Double.MAX_VALUE;
	
	for(TargetNode t1: p1.ap.values())
	{
		for(TargetNode t2: p2.ap.values())
		{
			
			TargetNode n1 =  SecurityGameContraction.getTargetNode(t1.getTargetid(), tmpgraph);
			TargetNode n2 = SecurityGameContraction.getTargetNode(t2.getTargetid(), tmpgraph);
			
			if(n1.getNeighbors().contains(n2) && (min>n1.getDistance(n2)))
			{
				min = t1.getDistance(t2);
			}
		}
	}
	
	
	return min;
}


private static double[] dOWithAttackCluster(int[][] gamedata,
		int nTargets, int nRes, double[][] density, double
		dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int RADIUS, HashMap<Integer, Integer> clusterhistogram, int slavelimit, int pathlimit) throws Exception {


	

	/*targets.clear();
	SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
	
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

	//printtargets(targets);

	/**
	 * 1. sort the targets
	 */
	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	SecurityGameContraction.printSortedTargets(targetssorted);
	
	
	HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>();
	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	

	HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
	HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
	
	
	
	int apspmat[][] =  buildAPSP(targets, nTargets, apspmap, apspmapback, apsp);
	
	

	ArrayList<Integer> currenttargets = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
	//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
	/*currenttargets.add(targetssorted[0][0]);
	currenttargets.add(targetssorted[1][0]);*/



	int currentPlace = currenttargets.size()-1;


	ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

	ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
	int attackedtarget=-1;
	ArrayList<Integer> attackhistory = new ArrayList<Integer>();
	int[][] p;
	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	HashSet jSet=new HashSet();
	ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
	List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
	double[] probdistribution= new double[jset.size()];
	double attackeru= -999;
	double attackerv = -999;

	


	long contractiontime=0;
	long solvingtime=0;
	long revmaptime=0;
	int targetsize=0;
	double slavetime = 0;
	int [][] origpmat = new int[nTargets][];
	int totalslaveiter = 0;
	long clusteringtime = 0;






	boolean canaddpath = true;
	boolean clusteringactivated = false;
	int masteritr=0;

	
	
	HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();
	
	
	
	// cluster history which will be maintained always
	HashMap<Integer, ArrayList<Integer>> attackclustershisotry = new HashMap<Integer, ArrayList<Integer>>();
	
	// list of targets which has been attacked
	ArrayList<Integer> clusteredtargets = new ArrayList<Integer>();
	
	// current attacked targets in defender oracle
	ArrayList<Integer> currentattackedtargets = new ArrayList<Integer>();
	
	// first target is the center of a cluster
	HashMap<Integer,Integer> clustercenters = new HashMap<Integer,Integer>();
	
	
	//cluster access points 
	HashMap<Integer, int[]> clusterap = new HashMap<Integer, int[]>();
	
	
	// current attacked Super targets in defender oracle
	ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();
	
	
	
	
	while(true)
	{
		
		
		
		
		System.out.println("\n clusteringactivated "+clusteringactivated +" Outer loop...Master");

		pathseq = new ArrayList<ArrayList<Integer>>();

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current place : "+ currentPlace);

		System.out.print("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}
		
		p = new int[targets.size()][]; // p matrix



		Date start = new Date();
		long l1 = start.getTime();

		
		tmpgraph.clear();
		domindatednodes.clear();
		HashMap<Integer, TargetNode> tmpgraphmaps = new HashMap<Integer, TargetNode>();
		
		generateGraph(tmpgraph, currentPlace, domindatednodes, targetssorted, dmax, targets, tmpgraphmaps);
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps);
		
		
		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		
		map.clear();
		mapback.clear();

		if(!clusteringactivated)
		{
			pathseq = generatePaths(tmpgraph, targets, dmax, currenttargets, currentPlace, pathseq, gamedata, mapback, map, nRes);
			
		}
		
		
		
		if(clusteringactivated)
		{
			System.out.println("need clustering "+clusteringactivated+", graph size "+ tmpgraph.size());
			//  build a graph with super targets
			// build an attack cluster except the dominated targets
			//sts.clear();
			dstravel.clear();
			stpaths.clear();
			
			
			
			start = new Date();
			l1 = start.getTime();

			


			
			sts = constructSuperTargets(tmpgraphmaps, attackhistory, apsp, apspmat, apspmap,apspmapback, dstravel, stpaths, (int) dmax, attackclustershisotry,
					clusteredtargets, clustercenters, currentattackedtargets, domindatednodes, RADIUS, tmpgraph, clusterap);
			
			

			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			clusteringtime += diff;
			

			//printSuperTargets(sts, stpaths, dstravel);
			preparePaths(dstravel, stpaths, sts);
			assignSTValues(sts, tmpgraphmaps);
			
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel);
			
			if(pathseq.size()==0)
			{
				System.out.println("No path seq................");
				//throw new Exception("No path seq");
			}
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths before removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths after removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: sts.values())
			{

				map.put(st.stid, icount);
				//System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			
			
			
			
			
			
		}

		

		int itr=0;
		currentattackedtargets.clear();
		currentattackedsupertargets.clear();
		while(true)
		{
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", Entered inner loop...slave");
			
			itr++;
			totalslaveiter++;

			

			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter.ohhh"+ iter);
			}*/


			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			
			
			if(!clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr);
					//SecurityGameContraction.printJointSchedule(jset);

					p = SecurityGameContraction.makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = SecurityGameContraction.findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = SecurityGameContraction.makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, gamedata);

					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" attack target after rev map "+ attackedtarget);
					
					
					//keep track of current attacked targets
					if(!currentattackedtargets.contains(attackedtarget))
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", "+ attackedtarget +" is added to the current attack list");
						currentattackedtargets.add(attackedtarget);
					}

					if(!attackhistory.contains(attackedtarget))
					{
						attackhistory.add(attackedtarget);
					}


					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
				{
					System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

					break;
				}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
				{
					System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
					break;
				}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);


					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;


					/**test
					 * 
					 */
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPathsDO(newpathseq, p, probdistribution, map);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" newpathseq size after purify : "+newpathseq.size());

					
					
					
					
					//////////////////FOR TESTING//////////////////
					clusteringactivated = true;
					if(clusteringactivated)
						break;
					/////////////////////////////////////
					
					

				} // end if else
			}// end if for clustering activated
			else if(clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: sts.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
					System.out.println("pathseq 0, iter"+ iter);
					if(currenttargets.size()==nTargets)
					{
						attackeru = mAxpayoff;
						attackerv = mAxpayoff;
						canaddpath = false;
						break;
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, sts, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, sts, tmpgraphmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, sts, map, mapback);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target before rev map "+ attackedtarget);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack ST before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, sts, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmpgraphmaps, sts, stpaths);
					
					
					
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmpgraphmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmpgraphmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target after rev map"+ attackedtarget);
					
					
					if(!currentattackedtargets.contains(attackedtarget))
					{
						currentattackedtargets.add(attackedtarget);
					}
					
					
					
					// find the attacked cluster
					// add it to the attackedsupertarget
					
					int hiddenattackedclsuter = findAttackedCluster(sts,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}
					
					
					
					
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +", inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();
					
					//printSuperTargets(sts);
					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmpgraphmaps, 
							dmax, sts.size(), 0, nRes, attackerstrategy, sts, dstravel, pathlimit);
					
					
					
					//////////// ADD Attacked clusters in the path////////
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);
					
					//////////////////////////////////////
					
					

					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=slavelimit))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" newpathseq: ");
						//SecurityGameContraction.printPaths(newpathseq);

						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Old path seq size "+ pathseq.size());

						//int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						//System.out.println("\n masteritr "+masteritr+" slaveitr "+itr +" new paths added by slave *************, attacked target "+ attackedtarget);

						//pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						//int newsize = pathseq.size();
						//System.out.println("haa ");


						/*if((oldsize==newsize) || (itr>=10))
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" Slave can't add any new path ############### or iteration>10");
							//printSuperTargets(sts);
						//	SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
			}
			System.out.println("iter"+ itr);
			
		} // inner while loop 




		// add all targets all targets with utility >= U(a')


		if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
			break;
		}

		


		double ulimit = SecurityGameContraction.getTargetNode(attackedtarget, targets).attackerreward;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


		int addcount=0;

		for(int k=currentPlace+1; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				addcount++;
				
				currenttargets.add(targetssorted[k][0]);
				System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
				if(addcount>=5)
				{
					break;
				}
			}
		}

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " addcount : "+ addcount);

		currentPlace = currenttargets.size()-1;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " currentplace  : "+ currentPlace);

		if(addcount<5 || addcount==0)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding more ");

			int prevcur = currentPlace;
			currentPlace += 5-addcount;

			//System.out.println("currentplace  : "+ currentPlace);
			if(currentPlace>targetssorted.length-1)
			{
				currentPlace = targetssorted.length-1;
			}
			//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{

				System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding target  "+ targetssorted[k][0]);
				currenttargets.add(targetssorted[k][0]);
			}
		}






		/*int prevcur = currentPlace;
		currentPlace += 3;
		if(currentPlace>targetssorted.length)
		{
			currentPlace = targetssorted.length;
		}
		System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

		for(int k= prevcur+1; k<=currentPlace; k++ )
		{
			currenttargets.add(targetssorted[k][0]);
		}


		break;*/
		masteritr++;


	} // outer while loop

	System.out.println("Final target list size : "+ currenttargets.size());

	for(int i=0; i<currenttargets.size(); i++)
	{
		System.out.print(currenttargets.get(i)+",");
	}

	//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
	double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);

	System.out.println("def exp : "+ defpayoff);
	
	
	//printST(sts, nTargets, iter);
	
	//printClusteredNodes(sts,nTargets,iter);
	
	//printClusterDists(sts,nTargets,iter);
	
	
	/*for(SuperTarget st: sts.values())
	{
		int value = st.nodes.size();
		if(clusterhistogram.containsKey(st.stid))
		{
			value += clusterhistogram.get(st.stid);
			
		}
		clusterhistogram.put(st.stid, value);
	}
	*/
	
	//ButtonGrid grid = new ButtonGrid(targetmaps, sts, "DOWithClus");
	//grid.drawPayoffGrid(nrow, ncol);
	//grid.drawCluster(nrow, ncol);


	//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

	double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime, totalslaveiter, clusteringtime};
	return res;
}



private static double[] dOWithAttackCluster2(int[][] gamedata,
		int nTargets, int nRes, double[][] density, double
		dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int RADIUS, HashMap<Integer, Integer> clusterhistogram, int slavelimit, int pathlimit) throws Exception {


	

	/*targets.clear();
	SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
	
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

	//printtargets(targets);

	/**
	 * 1. sort the targets
	 */
	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	SecurityGameContraction.printSortedTargets(targetssorted);
	
	
	HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>();
	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	

	HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
	HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
	
	
	
	int apspmat[][] =  buildAPSP(targets, nTargets, apspmap, apspmapback, apsp);
	
	

	ArrayList<Integer> currenttargets = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
	//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
	/*currenttargets.add(targetssorted[0][0]);
	currenttargets.add(targetssorted[1][0]);*/



	int currentPlace = currenttargets.size()-1;


	ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

	ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
	int attackedtarget=-1;
	ArrayList<Integer> attackhistory = new ArrayList<Integer>();
	int[][] p;
	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	HashSet jSet=new HashSet();
	ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
	List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
	double[] probdistribution= new double[jset.size()];
	double attackeru= -999;
	double attackerv = -999;

	


	long contractiontime=0;
	long solvingtime=0;
	long revmaptime=0;
	int targetsize=0;
	double slavetime = 0;
	int [][] origpmat = new int[nTargets][];
	int totalslaveiter = 0;
	long clusteringtime = 0;






	boolean canaddpath = true;
	boolean clusteringactivated = false;
	int masteritr=0;

	
	
	HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();
	
	
	
	// cluster history which will be maintained always
	HashMap<Integer, ArrayList<Integer>> attackclustershisotry = new HashMap<Integer, ArrayList<Integer>>();
	
	// list of targets which has been attacked
	ArrayList<Integer> clusteredtargets = new ArrayList<Integer>();
	
	// current attacked targets in defender oracle
	ArrayList<Integer> currentattackedtargets = new ArrayList<Integer>();
	
	// first target is the center of a cluster
	HashMap<Integer,Integer> clustercenters = new HashMap<Integer,Integer>();
	
	
	//cluster access points 
	HashMap<Integer, int[]> clusterap = new HashMap<Integer, int[]>();
	
	
	// current attacked Super targets in defender oracle
	ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();
	
	
	
	
	while(true)
	{
		
		
		
		
		System.out.println("\n clusteringactivated "+clusteringactivated +" Outer loop...Master");

		pathseq = new ArrayList<ArrayList<Integer>>();

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current place : "+ currentPlace);

		System.out.print("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}
		
		p = new int[targets.size()][]; // p matrix



		Date start = new Date();
		long l1 = start.getTime();

		
		tmpgraph.clear();
		domindatednodes.clear();
		HashMap<Integer, TargetNode> tmpgraphmaps = new HashMap<Integer, TargetNode>();
		
		generateGraph(tmpgraph, currentPlace, domindatednodes, targetssorted, dmax, targets, tmpgraphmaps);
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps);
		
		
		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		
		map.clear();
		mapback.clear();

		if(!clusteringactivated)
		{
			pathseq = generatePaths(tmpgraph, targets, dmax, currenttargets, currentPlace, pathseq, gamedata, mapback, map, nRes);
			
		}
		
		
		
		if(clusteringactivated)
		{
			System.out.println("need clustering "+clusteringactivated+", graph size "+ tmpgraph.size());
			//  build a graph with super targets
			// build an attack cluster except the dominated targets
			//sts.clear();
			dstravel.clear();
			stpaths.clear();
			
			
			
			start = new Date();
			l1 = start.getTime();

			


			
			sts = constructSuperTargets(tmpgraphmaps, attackhistory, apsp, apspmat, apspmap,apspmapback, dstravel, stpaths, (int) dmax, attackclustershisotry,
					clusteredtargets, clustercenters, currentattackedtargets, domindatednodes, RADIUS, tmpgraph, clusterap);
			
			

			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			clusteringtime += diff;
			

			//printSuperTargets(sts, stpaths, dstravel);
			preparePaths(dstravel, stpaths, sts);
			assignSTValues(sts, tmpgraphmaps);
			
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel);
			
			if(pathseq.size()==0)
			{
				System.out.println("No path seq................");
				//throw new Exception("No path seq");
			}
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths before removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths after removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: sts.values())
			{

				map.put(st.stid, icount);
				//System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			
			
			
			
			
			
		}

		

		int itr=0;
		currentattackedtargets.clear();
		currentattackedsupertargets.clear();
		while(true)
		{
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", Entered inner loop...slave");
			
			itr++;
			totalslaveiter++;

			

			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter.ohhh"+ iter);
			}*/


			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			
			
			if(!clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr);
					//SecurityGameContraction.printJointSchedule(jset);

					p = SecurityGameContraction.makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = SecurityGameContraction.findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = SecurityGameContraction.makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, gamedata);

					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" attack target after rev map "+ attackedtarget);
					
					
					//keep track of current attacked targets
					if(!currentattackedtargets.contains(attackedtarget))
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", "+ attackedtarget +" is added to the current attack list");
						currentattackedtargets.add(attackedtarget);
					}

					if(!attackhistory.contains(attackedtarget))
					{
						attackhistory.add(attackedtarget);
					}


					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
				{
					System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

					break;
				}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
				{
					System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
					break;
				}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);


					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;


					/**test
					 * 
					 */
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPathsDO(newpathseq, p, probdistribution, map);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" newpathseq size after purify : "+newpathseq.size());

					
					
					
					
					
					
					


					if(!clusteringactivated && (newpathseq.size()==0) /*|| (itr>=10)*/)
					{
						canaddpath = false;
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" Slave can't add any new path ###############");
						break;
					}
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" New path seq ");

				//	SecurityGameContraction.printPaths(newpathseq);








					//test
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);


					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);

					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+" new paths added by slave *************, attacked target "+ attackedtarget);

					/*pathseq = removeDuplicatePathSimple(pathseq);
				System.out.println("New path seq size "+ pathseq.size());
				//printPaths(pathseq);
				int newsize = pathseq.size();
				System.out.println("haa ");


				if((oldsize==newsize) || (itr>=10))
				{
					canaddpath = false;
					System.out.println("Slave can't add any new path ###############");
					break;
				}*/
					
					
					
					if(!clusteringactivated && newpathseq.size()>0 && itr>=slavelimit)
					{
						
						clusteringactivated = true;
						System.out.println("\n clusteringactivated "+clusteringactivated +"  masteritr "+masteritr+ ", slaveitr "+itr+" Need clustering "+ clusteringactivated);
						break;
					}

					System.out.println("\n clusteringactivated "+clusteringactivated +"  masteritr "+masteritr+ ", slaveitr "+itr+" New full path seq after adding new paths");

					//SecurityGameContraction.printPaths(pathseq);


				} // end if else
			}// end if for clustering activated
			else if(clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: sts.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
					System.out.println("pathseq 0, iter"+ iter);
					if(currenttargets.size()==nTargets)
					{
						attackeru = mAxpayoff;
						attackerv = mAxpayoff;
						canaddpath = false;
						break;
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, sts, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, sts, tmpgraphmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, sts, map, mapback);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target before rev map "+ attackedtarget);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack ST before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, sts, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmpgraphmaps, sts, stpaths);
					
					
					
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmpgraphmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmpgraphmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target after rev map"+ attackedtarget);
					
					
					if(!currentattackedtargets.contains(attackedtarget))
					{
						currentattackedtargets.add(attackedtarget);
					}
					
					
					
					// find the attacked cluster
					// add it to the attackedsupertarget
					
					int hiddenattackedclsuter = findAttackedCluster(sts,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}
					
					
					
					
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +", inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();
					
					//printSuperTargets(sts);
					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmpgraphmaps, 
							dmax, sts.size(), 0, nRes, attackerstrategy, sts, dstravel, pathlimit);
					
					
					
					//////////// ADD Attacked clusters in the path////////
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);
					
					//////////////////////////////////////
					
					

					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=slavelimit))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" newpathseq: ");
						//SecurityGameContraction.printPaths(newpathseq);

						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Old path seq size "+ pathseq.size());

						//int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						//System.out.println("\n masteritr "+masteritr+" slaveitr "+itr +" new paths added by slave *************, attacked target "+ attackedtarget);

						//pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						//int newsize = pathseq.size();
						//System.out.println("haa ");


						/*if((oldsize==newsize) || (itr>=10))
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" Slave can't add any new path ############### or iteration>10");
							//printSuperTargets(sts);
						//	SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
			}
			System.out.println("iter"+ itr);
			
		} // inner while loop 




		// add all targets all targets with utility >= U(a')


		if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
			break;
		}

		


		double ulimit = SecurityGameContraction.getTargetNode(attackedtarget, targets).attackerreward;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


		int addcount=0;

		for(int k=currentPlace+1; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				addcount++;
				
				currenttargets.add(targetssorted[k][0]);
				System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
				if(addcount>=5)
				{
					break;
				}
			}
		}

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " addcount : "+ addcount);

		currentPlace = currenttargets.size()-1;

		System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " currentplace  : "+ currentPlace);

		if(addcount<5 || addcount==0)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding more ");

			int prevcur = currentPlace;
			currentPlace += 5-addcount;

			//System.out.println("currentplace  : "+ currentPlace);
			if(currentPlace>targetssorted.length-1)
			{
				currentPlace = targetssorted.length-1;
			}
			//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{

				System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding target  "+ targetssorted[k][0]);
				currenttargets.add(targetssorted[k][0]);
			}
		}






		/*int prevcur = currentPlace;
		currentPlace += 3;
		if(currentPlace>targetssorted.length)
		{
			currentPlace = targetssorted.length;
		}
		System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

		for(int k= prevcur+1; k<=currentPlace; k++ )
		{
			currenttargets.add(targetssorted[k][0]);
		}


		break;*/
		masteritr++;


	} // outer while loop

	System.out.println("Final target list size : "+ currenttargets.size());

	for(int i=0; i<currenttargets.size(); i++)
	{
		System.out.print(currenttargets.get(i)+",");
	}

	//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
	double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);

	System.out.println("def exp : "+ defpayoff);
	
	
	//printST(sts, nTargets, iter);
	
	//printClusteredNodes(sts,nTargets,iter);
	
	//printClusterDists(sts,nTargets,iter);
	
	
	/*for(SuperTarget st: sts.values())
	{
		int value = st.nodes.size();
		if(clusterhistogram.containsKey(st.stid))
		{
			value += clusterhistogram.get(st.stid);
			
		}
		clusterhistogram.put(st.stid, value);
	}
	*/
	
	//ButtonGrid grid = new ButtonGrid(targetmaps, sts, "DOWithClus");
	//grid.drawPayoffGrid(nrow, ncol);
	//grid.drawCluster(nrow, ncol);


	//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

	double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime, totalslaveiter, clusteringtime};
	return res;
}



private static double[] dOWithAttackCluster3(int[][] gamedata,
		int nTargets, int nRes, double[][] density, double
		dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int RADIUS, HashMap<Integer, Integer> clusterhistogram, int slavelimit, int pathlimit) throws Exception {


	

	/*targets.clear();
	SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
	
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

	//printtargets(targets);

	/**
	 * 1. sort the targets
	 */
	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	SecurityGameContraction.printSortedTargets(targetssorted);
	
	
	HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>();
	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	

	HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
	HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
	
	
	
	int apspmat[][] =  buildAPSP(targets, nTargets, apspmap, apspmapback, apsp);
	
	

	ArrayList<Integer> currenttargets = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
	//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
	/*currenttargets.add(targetssorted[0][0]);
	currenttargets.add(targetssorted[1][0]);*/



	int currentPlace = currenttargets.size()-1;


	ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

	ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
	int attackedtarget=-1;
	ArrayList<Integer> attackhistory = new ArrayList<Integer>();
	int[][] p;
	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	HashSet jSet=new HashSet();
	ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
	List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
	double[] probdistribution= new double[jset.size()];
	double attackeru= -999;
	double attackerv = -999;

	


	long contractiontime=0;
	long solvingtime=0;
	long revmaptime=0;
	int targetsize=0;
	double slavetime = 0;
	int [][] origpmat = new int[nTargets][];
	int totalslaveiter = 0;
	long clusteringtime = 0;






	boolean canaddpath = true;
	boolean clusteringactivated = false;
	int masteritr=0;

	
	
	HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();
	
	
	
	// cluster history which will be maintained always
	HashMap<Integer, ArrayList<Integer>> attackclustershisotry = new HashMap<Integer, ArrayList<Integer>>();
	
	// list of targets which has been attacked
	ArrayList<Integer> clusteredtargets = new ArrayList<Integer>();
	
	// current attacked targets in defender oracle
	ArrayList<Integer> currentattackedtargets = new ArrayList<Integer>();
	
	// first target is the center of a cluster
	HashMap<Integer,Integer> clustercenters = new HashMap<Integer,Integer>();
	
	
	//cluster access points 
	HashMap<Integer, int[]> clusterap = new HashMap<Integer, int[]>();
	
	
	// current attacked Super targets in defender oracle
	ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();
	
	clusteringactivated = true;
	
	
	ArrayList<Double[]> masterslaveres = new ArrayList<Double[]>();
	
	
	
	while(true)
	{
		
		
		
		
		System.out.println("\n clusteringactivated "+clusteringactivated +" Outer loop...Master");

		pathseq = new ArrayList<ArrayList<Integer>>();

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current place : "+ currentPlace);

		//System.out.print("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current target list : ");

		/*for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}*/
		
		p = new int[targets.size()][]; // p matrix



		Date start = new Date();
		long l1 = start.getTime();

		
		tmpgraph.clear();
		domindatednodes.clear();
		HashMap<Integer, TargetNode> tmpgraphmaps = new HashMap<Integer, TargetNode>();
		
		generateGraph(tmpgraph, currentPlace, domindatednodes, targetssorted, dmax, targets, tmpgraphmaps);
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps);
		
		
		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		
		map.clear();
		mapback.clear();

		if(clusteringactivated)
		{
			System.out.println("need clustering "+clusteringactivated+", graph size "+ tmpgraph.size());
			//  build a graph with super targets
			// build an attack cluster except the dominated targets
			//sts.clear();
			dstravel.clear();
			stpaths.clear();
			
			
			
			start = new Date();
			l1 = start.getTime();

			


			
			sts = constructSuperTargets(tmpgraphmaps, attackhistory, apsp, apspmat, apspmap,apspmapback, dstravel, stpaths, (int) dmax, attackclustershisotry,
					clusteredtargets, clustercenters, currentattackedtargets, domindatednodes, RADIUS, tmpgraph, clusterap);
			
			

			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			clusteringtime += diff;
			

			printSuperTargets(sts, stpaths, dstravel);
			preparePaths(dstravel, stpaths, sts);
			assignSTValues(sts, tmpgraphmaps);
			
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel);
			
			if(pathseq.size()==0)
			{
				System.out.println("No path seq................");
				//throw new Exception("No path seq");
			}
			
			//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths before removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths after removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: sts.values())
			{

				map.put(st.stid, icount);
				//System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			
			
			
			
			
			
		}

		

		int itr=0;
		currentattackedtargets.clear();
		currentattackedsupertargets.clear();
		while(true)
		{
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", Entered inner loop...slave");
			
			itr++;
			totalslaveiter++;

			

			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter.ohhh"+ iter);
			}*/


			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			
			
			if(clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: sts.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
					System.out.println("pathseq 0, iter"+ iter);
					if(currenttargets.size()==nTargets)
					{
						attackeru = mAxpayoff;
						attackerv = mAxpayoff;
						canaddpath = false;
						break;
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, sts, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, sts, null, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, sts, map, mapback);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target before rev map "+ attackedtarget);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack ST before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, sts, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmpgraphmaps, sts, stpaths);
					
					
					
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmpgraphmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmpgraphmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target after rev map"+ attackedtarget);
					
					
					if(!currentattackedtargets.contains(attackedtarget))
					{
						currentattackedtargets.add(attackedtarget);
					}
					
					
					
					// find the attacked cluster
					// add it to the attackedsupertarget
					
					int hiddenattackedclsuter = findAttackedCluster(sts,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}
					
					
					
					
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					/*if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}
*/
					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +", inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();
					
					//printSuperTargets(sts);
					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmpgraphmaps, dmax, sts.size(), 0, nRes, attackerstrategy, sts, dstravel, pathlimit);
					
					//ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();
					
					//////////// ADD Attacked clusters in the path////////
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);
					
					//////////////////////////////////////
					
					

					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					//System.out.println("newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
					//System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=slavelimit))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						//System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" newpathseq: ");
						//SecurityGameContraction.printPaths(newpathseq);

						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Old path seq size "+ pathseq.size());

						//int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						//System.out.println("\n masteritr "+masteritr+" slaveitr "+itr +" new paths added by slave *************, attacked target "+ attackedtarget);

						//pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						//int newsize = pathseq.size();
						//System.out.println("haa ");


						/*if((oldsize==newsize) || (itr>=10))
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" Slave can't add any new path ############### or iteration>10");
							//printSuperTargets(sts);
						//	SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
			}
			System.out.println("iter"+ itr);
			
			Double[] key = {(double)masteritr ,(double)itr-1, attackerv};
			masterslaveres.add(key);
			
		} // inner while loop 




		// add all targets all targets with utility >= U(a')


		if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
			break;
		}

		


		double ulimit = SecurityGameContraction.getTargetNode(attackedtarget, targets).attackerreward;

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


		int addcount=0;

		for(int k=currentPlace+1; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				addcount++;
				
				currenttargets.add(targetssorted[k][0]);
				//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
				if(addcount>=5)
				{
					break;
				}
			}
		}

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " addcount : "+ addcount);

		currentPlace = currenttargets.size()-1;

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " currentplace  : "+ currentPlace);

		if(addcount<5 || addcount==0)
		{
			//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding more ");

			int prevcur = currentPlace;
			currentPlace += 5-addcount;

			//System.out.println("currentplace  : "+ currentPlace);
			if(currentPlace>targetssorted.length-1)
			{
				currentPlace = targetssorted.length-1;
			}
			//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{

				//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding target  "+ targetssorted[k][0]);
				currenttargets.add(targetssorted[k][0]);
			}
		}






		/*int prevcur = currentPlace;
		currentPlace += 3;
		if(currentPlace>targetssorted.length)
		{
			currentPlace = targetssorted.length;
		}
		System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

		for(int k= prevcur+1; k<=currentPlace; k++ )
		{
			currenttargets.add(targetssorted[k][0]);
		}


		break;*/
		masteritr++;


	} // outer while loop

	System.out.println("Final target list size : "+ currenttargets.size());

	/*for(int i=0; i<currenttargets.size(); i++)
	{
		System.out.print(currenttargets.get(i)+",");
	}*/

	//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
	double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);

	System.out.println("def exp : "+ defpayoff);
	
	
	//printST(sts, nTargets, iter);
	
	//printClusteredNodes(sts,nTargets,iter);
	
	//printClusterDists(sts,nTargets,iter);
	
	/*
	for(SuperTarget st: sts.values())
	{
		int value = st.nodes.size();
		if(clusterhistogram.containsKey(st.stid))
		{
			value += clusterhistogram.get(st.stid);
			
		}
		clusterhistogram.put(st.stid, value);
	}*/
	
	
	//ButtonGrid grid = new ButtonGrid(targetmaps, sts, "DOWithClus");
	//grid.drawPayoffGrid(nrow, ncol);
	//grid.drawCluster(nrow, ncol);


	//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
	
	//printSuperTargets(sts);
	
	//verifySolution(jset, pathseq, probdistribution, nTargets, dmax, sts, dstravel, tmpgraph);
	
	//writeMasterSlaveRes(masterslaveres);
	

	double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime, totalslaveiter, clusteringtime};
	return res;
}




/**
 * attack clustering without contraction
 * @param gamedata
 * @param nTargets
 * @param nRes
 * @param density
 * @param dmax
 * @param iter
 * @param nrow
 * @param ncol
 * @param targets
 * @param targetmaps
 * @param RADIUS
 * @param clusterhistogram
 * @param slavelimit
 * @param pathlimit
 * @return
 * @throws Exception
 */
private static double[] dOWithAttackCluster4(int[][] gamedata,
		int nTargets, int nRes, double[][] density, double
		dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int RADIUS, HashMap<Integer, Integer> clusterhistogram, int slavelimit, int pathlimit) throws Exception {


	

	/*targets.clear();
	SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
	
	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	
	assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

	//printtargets(targets);

	/**
	 * 1. sort the targets
	 */
	int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
	SecurityGameContraction.printSortedTargets(targetssorted);
	
	
	HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>();
	AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
	

	HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
	HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
	
	
	
	int apspmat[][] =  buildAPSP(targets, nTargets, apspmap, apspmapback, apsp);
	
	

	ArrayList<Integer> currenttargets = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
	//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
	/*currenttargets.add(targetssorted[0][0]);
	currenttargets.add(targetssorted[1][0]);*/



	int currentPlace = currenttargets.size()-1;


	ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

	ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
	int attackedtarget=-1;
	ArrayList<Integer> attackhistory = new ArrayList<Integer>();
	int[][] p;
	HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
	HashSet jSet=new HashSet();
	ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
	List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
	double[] probdistribution= new double[jset.size()];
	double attackeru= -999;
	double attackerv = -999;

	


	long contractiontime=0;
	long solvingtime=0;
	long revmaptime=0;
	int targetsize=0;
	double slavetime = 0;
	int [][] origpmat = new int[nTargets][];
	int totalslaveiter = 0;
	long clusteringtime = 0;






	boolean canaddpath = true;
	boolean clusteringactivated = false;
	int masteritr=0;

	
	
	HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();
	
	
	
	// cluster history which will be maintained always
	HashMap<Integer, ArrayList<Integer>> attackclustershisotry = new HashMap<Integer, ArrayList<Integer>>();
	
	// list of targets which has been attacked
	ArrayList<Integer> clusteredtargets = new ArrayList<Integer>();
	
	// current attacked targets in defender oracle
	ArrayList<Integer> currentattackedtargets = new ArrayList<Integer>();
	
	// first target is the center of a cluster
	HashMap<Integer,Integer> clustercenters = new HashMap<Integer,Integer>();
	
	
	//cluster access points 
	HashMap<Integer, int[]> clusterap = new HashMap<Integer, int[]>();
	
	
	// current attacked Super targets in defender oracle
	ArrayList<Integer> currentattackedsupertargets = new ArrayList<Integer>();
	
	clusteringactivated = true;
	
	
	ArrayList<Double[]> masterslaveres = new ArrayList<Double[]>();
	
	
	
	while(true)
	{
		
		
		
		
		System.out.println("\n clusteringactivated "+clusteringactivated +" Outer loop...Master");

		pathseq = new ArrayList<ArrayList<Integer>>();

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current place : "+ currentPlace);

		//System.out.print("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Current target list : ");

		/*for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}*/
		
		p = new int[targets.size()][]; // p matrix



		Date start = new Date();
		long l1 = start.getTime();

		
		tmpgraph.clear();
		domindatednodes.clear();
		HashMap<Integer, TargetNode> tmpgraphmaps = new HashMap<Integer, TargetNode>();
		
		generateGraphWithOutContraction(tmpgraph, currentPlace, domindatednodes, targetssorted, dmax, targets, tmpgraphmaps);
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps, domindatednodes);
	
		
		
		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		
		map.clear();
		mapback.clear();

		if(clusteringactivated)
		{
			System.out.println("need clustering "+clusteringactivated+", graph size "+ tmpgraph.size());
			//  build a graph with super targets
			// build an attack cluster except the dominated targets
			//sts.clear();
			dstravel.clear();
			stpaths.clear();
			
			
			
			start = new Date();
			l1 = start.getTime();

			

			// dominated targets will be single cluster/ super target
			
			sts = constructSuperTargetsV3(tmpgraphmaps, attackhistory, apsp, apspmat, apspmap,apspmapback, dstravel, stpaths, (int) dmax, attackclustershisotry,
					clusteredtargets, clustercenters, currentattackedtargets, domindatednodes, RADIUS, tmpgraph, clusterap);
			
			

			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			clusteringtime += diff;
			

			printSuperTargets(sts, stpaths, dstravel);
			preparePaths(dstravel, stpaths, sts);
			assignSTValues(sts, tmpgraphmaps);
			
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel);
			
			if(pathseq.size()==0)
			{
				System.out.println("No path seq................");
				//throw new Exception("No path seq");
			}
			
			//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths before removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" Paths after removing duplicate : ");
			//SecurityGameContraction.printPaths(pathseq);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: sts.values())
			{

				map.put(st.stid, icount);
				//System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			
			
			
			
			
			
		}

		

		int itr=0;
		currentattackedtargets.clear();
		currentattackedsupertargets.clear();
		while(true)
		{
			
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ ", slaveitr "+itr+", Entered inner loop...slave");
			
			itr++;
			totalslaveiter++;

			

			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter.ohhh"+ iter);
			}*/


			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			
			
			if(clusteringactivated)
			{
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: sts.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}
					}
					System.out.println("pathseq 0, iter"+ iter);
					if(currenttargets.size()==nTargets)
					{
						attackeru = mAxpayoff;
						attackerv = mAxpayoff;
						canaddpath = false;
						break;
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, sts, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, sts, null, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, sts, map, mapback);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target before rev map "+ attackedtarget);
					attackedtarget = mapback.get(attackedtarget);
					
					int attackedclsuter = attackedtarget;
					
					if(!currentattackedsupertargets.contains(attackedclsuter))
					{
						currentattackedsupertargets.add(attackedclsuter);
					}
					
					
					
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack ST before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, sts, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmpgraphmaps, sts, stpaths);
					
					
					
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmpgraphmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmpgraphmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("\n clusteringactivated "+clusteringactivated +" master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" attack target after rev map"+ attackedtarget);
					
					
					if(!currentattackedtargets.contains(attackedtarget))
					{
						currentattackedtargets.add(attackedtarget);
					}
					
					
					
					// find the attacked cluster
					// add it to the attackedsupertarget
					
					int hiddenattackedclsuter = findAttackedCluster(sts,attackedtarget);
					
					// might attack target which are dominated
					
					if(hiddenattackedclsuter != -1)
					{
						hiddenattackedclsuter = mapback.get(hiddenattackedclsuter);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" hidden attack ST "+ hiddenattackedclsuter);
						if(!currentattackedsupertargets.contains(hiddenattackedclsuter))
						{
							currentattackedsupertargets.add(hiddenattackedclsuter);
						}
					}
					
					
					
					
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					/*if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}
*/
					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +", inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();
					
					//printSuperTargets(sts);
					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSTSlavePaths(tmpgraphmaps, dmax, sts.size(), 0, nRes, attackerstrategy, sts, dstravel, pathlimit);
					
					//ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();
					
					//////////// ADD Attacked clusters in the path////////
					
					SecurityGameContraction.addSuperTargetsAPSP(dmax, sts, tmpgraphmaps, nRes, dstravel, newpathseq, currentattackedsupertargets);
					
					//////////////////////////////////////
					
					

					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					//System.out.println("newpathseq size before purify : "+newpathseq.size());
					newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
					//System.out.println("newpathseq size after purify : "+newpathseq.size());
						
						
						if((newpathseq.size()==0) || (itr>=slavelimit))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}
						//System.out.println("New whole path seq ");
						
						
						

						//makeSlavePathSeq(newpathseq, goal);
						//removeDuplicatePathSimple(newpathseq);
						/*if(newpathseq.size()==0)
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Slave can't add any new path ###############");
							break;
						}*/
						//System.out.println("tcur: ");
						//printGreedyPath(currenttargets);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" newpathseq: ");
						//SecurityGameContraction.printPaths(newpathseq);

						System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" Old path seq size "+ pathseq.size());

						//int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						//System.out.println("\n masteritr "+masteritr+" slaveitr "+itr +" new paths added by slave *************, attacked target "+ attackedtarget);

						//pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+" slaveitr "+itr +" New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						//int newsize = pathseq.size();
						//System.out.println("haa ");


						/*if((oldsize==newsize) || (itr>=10))
						{
							canaddpath = false;
							System.out.println("\n clusteringactivated "+clusteringactivated +" Slave can't add any new path ############### or iteration>10");
							//printSuperTargets(sts);
						//	SecurityGameContraction.printPaths(pathseq);
							break;
						}*/

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
			}
			System.out.println("iter"+ itr);
			
			Double[] key = {(double)masteritr ,(double)itr-1, attackerv};
			masterslaveres.add(key);
			
		} // inner while loop 




		// add all targets all targets with utility >= U(a')


		if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
		{
			System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
			break;
		}

		


		double ulimit = SecurityGameContraction.getTargetNode(attackedtarget, targets).attackerreward;

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


		int addcount=0;

		for(int k=currentPlace+1; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				addcount++;
				
				currenttargets.add(targetssorted[k][0]);
				//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
				if(addcount>=5)
				{
					break;
				}
			}
		}

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " addcount : "+ addcount);

		currentPlace = currenttargets.size()-1;

		//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " currentplace  : "+ currentPlace);

		if(addcount<5 || addcount==0)
		{
			//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding more ");

			int prevcur = currentPlace;
			currentPlace += 5-addcount;

			//System.out.println("currentplace  : "+ currentPlace);
			if(currentPlace>targetssorted.length-1)
			{
				currentPlace = targetssorted.length-1;
			}
			//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{

				//System.out.println("\n clusteringactivated "+clusteringactivated +" masteritr "+masteritr+ " adding target  "+ targetssorted[k][0]);
				currenttargets.add(targetssorted[k][0]);
			}
		}






		/*int prevcur = currentPlace;
		currentPlace += 3;
		if(currentPlace>targetssorted.length)
		{
			currentPlace = targetssorted.length;
		}
		System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

		for(int k= prevcur+1; k<=currentPlace; k++ )
		{
			currenttargets.add(targetssorted[k][0]);
		}


		break;*/
		masteritr++;


	} // outer while loop

	System.out.println("Final target list size : "+ currenttargets.size());

	/*for(int i=0; i<currenttargets.size(); i++)
	{
		System.out.print(currenttargets.get(i)+",");
	}*/

	//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
	double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);

	System.out.println("def exp : "+ defpayoff);
	
	
	//printST(sts, nTargets, iter);
	
	//printClusteredNodes(sts,nTargets,iter);
	
	//printClusterDists(sts,nTargets,iter);
	
	/*
	for(SuperTarget st: sts.values())
	{
		int value = st.nodes.size();
		if(clusterhistogram.containsKey(st.stid))
		{
			value += clusterhistogram.get(st.stid);
			
		}
		clusterhistogram.put(st.stid, value);
	}*/
	
	
	//ButtonGrid grid = new ButtonGrid(targetmaps, sts, "DOWithClus");
	//grid.drawPayoffGrid(nrow, ncol);
	//grid.drawCluster(nrow, ncol);


	//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
	
	//printSuperTargets(sts);
	
	//verifySolution(jset, pathseq, probdistribution, nTargets, dmax, sts, dstravel, tmpgraph);
	
	//writeMasterSlaveRes(masterslaveres);
	

	double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime, totalslaveiter, clusteringtime};
	return res;
}





private static void printNodesWithNeighborsAndPath(HashMap<Integer, TargetNode> tmpgraphmaps,
		ArrayList<TargetNode> domindatednodes) {
	
	for(TargetNode node : tmpgraphmaps.values())
	{
		
		//Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
		if(!domindatednodes.contains(node))
		{
			
			System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");
			for(TargetNode neighbor: node.getNeighbors())
			{
				
				if(node.getTargetid()==224 && neighbor.getTargetid()==272)
				{
					System.out.println();
				}
				
				System.out.println("---Neighbor : "+ neighbor.getTargetid());
				//Logger.logit("---Neighbor : "+ neighbor.getTargetid()+"\n");
				/**
				 * print path
				 */
				ArrayList<TargetNode> path = node.getPath(neighbor);
				System.out.print("Path : "+ node.getTargetid()+ " --> ");
				//Logger.logit("Path : "+ node.getTargetid()+ " --> ");
				for(TargetNode pathnode : path)
				{
					System.out.print(pathnode.getTargetid()+" --> ");
					Logger.logit(pathnode.getTargetid()+" --> ");
				}

				System.out.print(neighbor.getTargetid()+ "\n");
				System.out.println("Distance : " + node.getDistance(neighbor));
				//Logger.logit(neighbor.getTargetid()+ "\n");
				//Logger.logit("Distance : " + node.getDistance(neighbor));

				System.out.print("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");
				//Logger.logit("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");

			}
		}
	}
	
}


private static void writeMasterSlaveRes(ArrayList<Double[]> masterslaveres) {
	
	
	
	try
	{
		PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result/master-slave-result.csv"),true));
		//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
		
		for(Double[] key: masterslaveres)
		{
			//String k[] = key.split(",");
			//double value = masterslaveres.get(key);
		
			pw.append(key[0]+","+ key[1]+","+key[2]+"\n");
		}
		pw.close();

	}
	catch(Exception e)
	{

	}
	
}
	
	private static void printClusterDists(HashMap<Integer, SuperTarget> sts, int nTargets, int iter) {
	
		
		
			try
			{
				
				File f = new File("result\\clusdist"+nTargets+"-"+iter+".csv");
				 
				 if(f.exists())
				 {
					 f.delete();
					 f.createNewFile();
				 }
				
				
				PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result\\clusdist"+nTargets+"-"+iter+".csv"),true));
				
				for(SuperTarget st: sts.values())
				{

					

					if(st.nodes.size()>1)
					{
						pw.append("Cluster "+st.stid+"\n");
						for(TargetNode t: st.nodes.values())
						{
							for(TargetNode t2: st.nodes.values())
							{
								if(t.getNeighbors().contains(t2))
								{
									pw.append(t.getTargetid()+"->"+t2.getTargetid()+ "="+t.getDistance(t2)+"  Path : ");
									
									for(TargetNode pnode: t.getPath(t2))
									{
										pw.append(pnode.getTargetid()+"->");
									}
									pw.append(",");
								}
								
								
								//pw.append(t.getTargetid()+","+t.defenderreward+","+t.defenderpenalty+","+t.attackerreward+","+t.attackerpenalty+"\n");
							}
							pw.append("\n");

						}
						pw.append("\n");
					}

				}
				
				//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
				//pw.append(expno+","+nTargets+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," +slavetime+","+ totaltime+"\n");
				pw.close();

			}
			catch(Exception e)
			{

			}
			
		
	
	
}


	private static void printClusteredNodes(HashMap<Integer, SuperTarget> sts, int nTargets, int iter) {
	
		
			try
			{
				
				File f = new File("result\\clusnodefeatures"+nTargets+"-"+iter+".csv");
				 
				 if(f.exists())
				 {
					 f.delete();
					 f.createNewFile();
				 }
				
				
				PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result\\clusnodefeatures"+nTargets+"-"+iter+".csv"),true));
				
				pw.append("ID,DR,DP,AR,AP"+"\n");
				
				for(SuperTarget st: sts.values())
				{

					if(st.nodes.size()>1)
					{
						for(TargetNode t: st.nodes.values())
						{
							pw.append(t.getTargetid()+","+t.defenderreward+","+t.defenderpenalty+","+t.attackerreward+","+t.attackerpenalty+"\n");

						}
						pw.append("\n");
					}

				}
				
				//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
				//pw.append(expno+","+nTargets+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," +slavetime+","+ totaltime+"\n");
				pw.close();

			}
			catch(Exception e)
			{

			}
			
		
	
}


	private static void printST(HashMap<Integer, SuperTarget> sts, int nTargets, int iter) {
	// TODO Auto-generated method stub
		
		
		
			try
			{
				
				File f = new File("result\\clusters"+nTargets+"-"+iter+".csv");
				 
				 if(f.exists())
				 {
					 f.delete();
					 f.createNewFile();
				 }
				
				
				PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result\\clusters"+nTargets+"-"+iter+".csv"),true));
				
				for(SuperTarget st: sts.values())
				{
					
					pw.append("Cluster "+ st.stid + ": ,");
					for(TargetNode t: st.nodes.values())
					{
						pw.append(t.getTargetid()+",");

					}
					pw.append("\n");
				}
				
				
				
				
				//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
				//pw.append(expno+","+nTargets+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," +slavetime+","+ totaltime+"\n");
				pw.close();

			}
			catch(Exception e)
			{

			}
			
		
		
	
}


	private static int findAttackedCluster(HashMap<Integer, SuperTarget> sts, int attackedtarget) {
	
		
		for(SuperTarget st: sts.values())
		{
			for(TargetNode t: st.nodes.values())
			{
				if(t.getTargetid()  == attackedtarget)
				{
					return st.stid;
				}
			}
		}
		
		
		
	return -1;
}


	private static void preparePaths(HashMap<Integer, Double> dstravel, HashMap<Integer, ArrayList<Integer>> stpaths,
		HashMap<Integer, SuperTarget> sts) {
		
		
		ArrayList<Integer> notin = new ArrayList<Integer>();
		
		for(SuperTarget st: sts.values())
		{
			if(st.nodes.size()==1)
			{
				dstravel.put(st.stid, 0.0);
				
			}
			
		}
		
		
		int ind = 0;
		for(Integer t: dstravel.keySet())
		{
			if(!sts.keySet().contains(t))
			{
				notin.add(t);
			}
			ind++;
		}
		
		for(Integer x: notin)
		{
			dstravel.remove(x);
			stpaths.remove(x);
		}
	
}


	private static int[][] buildAPSP(ArrayList<TargetNode> targets, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback, AllPairShortestPath apsp) {
	

		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];
		
		
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		SecurityGameContraction.purifyAPSPMatrixZeroGT(apspmat, targets, nTargets, map, mapback);
		
		return apspmat;
	
}


	private static HashMap<Integer, SuperTarget> createSuperTargets(HashMap<Integer, TargetNode> tmpgraphmaps,
		int attackedtarget, ArrayList<Integer> attackhistory, AllPairShortestPath apsp, int[][] apspmat,
		HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback,
		HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax) {
	
		
		// for now build a cluster
		
		//for(Integer attackt: attackhistory)
		
		ArrayList<Integer> attackcluster = createAttackCluster(tmpgraphmaps, attackedtarget);
		
		int totalcluster = tmpgraphmaps.size() - attackcluster.size() + 1;
		
		
		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[totalcluster];
		
		for(int i=0; i<totalcluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}
		
		clusters[0].add(0);
		
		for(Integer n: attackcluster)
		{
		
			clusters[1].add(n);
		}
		
		int index = 2;
		for(Integer t: tmpgraphmaps.keySet())
		{
			if(t != 0 && !attackcluster.contains(t))
			{
				clusters[index++].add(t);
				
			}
		}
		
		
		
		
		printClusters(clusters);
		
		
		
		

		
		HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,tmpgraphmaps);
		
		printSuperTargets(sts);
		
		
		System.out.println("\n \n \n After choosing AP \n \n \n");
		
		
		chooseAP(sts, tmpgraphmaps, sts.get(1), apsp, apspmat, apspmap, apspmapback, dstravel, stpaths, dmax, null,null);
		
		
		return sts;
		
			
		
}
	
	
	
	private static HashMap<Integer, SuperTarget> constructSuperTargets(HashMap<Integer, TargetNode> tmpgraphmaps,
			 ArrayList<Integer> attackhistory, AllPairShortestPath apsp, int[][] apspmat,
			HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback,
			HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax,
			HashMap<Integer, ArrayList<Integer>> attackclustershisotry, ArrayList<Integer> clusteredtargets, 
			HashMap<Integer,Integer> clustercenters, ArrayList<Integer> newtargetstocluster, ArrayList<TargetNode> domindatednodes, 
			int RADIUS, ArrayList<TargetNode> tmpgraph, HashMap<Integer, int[]> clusterap) {
		
			
			// for now build a cluster
			
			//for(Integer attackt: attackhistory)
			
			//ArrayList<Integer> attackcluster = new ArrayList<Integer>();//createAttackCluster(tmpgraphmaps, attackedtarget);
			
			// base target will alws be alone, so no clustering using base target
		
		
		
		

		/*System.out.println("Current attack clusters ");


		printClusters(attackclustershisotry);

		System.out.println("Current attacked targets ");

		printTargets(newtargetstocluster);


		System.out.println("Clustered targets ");

		printTargets(clusteredtargets);
*/




		System.out.println("Updating attack clusters ");
		
		
		

		ArrayList<Integer> changedcluster = updateClusters(attackclustershisotry, clusteredtargets, clustercenters, newtargetstocluster, tmpgraphmaps, RADIUS, domindatednodes, tmpgraph);







		// create clusters

		int totalcluster = tmpgraphmaps.size() - clusteredtargets.size() + attackclustershisotry.size() - domindatednodes.size();


		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[totalcluster];

		for(int i=0; i<totalcluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}

		clusters[0].add(0);

		int clusindex = 1;
		for(ArrayList<Integer> attackcluster: attackclustershisotry.values())
		{

			for(Integer n: attackcluster)
			{

				clusters[clusindex].add(n);
			}
			clusindex++;
		}

		//int index = 2;
		for(Integer t: tmpgraphmaps.keySet())
		{
			TargetNode tnode = tmpgraphmaps.get(t);
			
			if(t != 0 && !clusteredtargets.contains(t) && (!domindatednodes.contains(tnode)))
			{
				clusters[clusindex++].add(t);

			}
		}

			
			
			
			//printClusters(clusters);
			
			
			
			

			
			HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,tmpgraphmaps);
			
			//printSuperTargets(sts);
			
			
			System.out.println("\n \n \n After choosing AP \n \n \n");
			
			
			for(int i = 0; i<sts.size(); i++)
			{
				SuperTarget st = sts.get(i);
				
				
				
				
				// if cluster is chnaged then compute AP
				// if not changed then use the previous AP
				
				try {
					chooseDegreeBasedAPV2(sts, tmpgraphmaps, sts.get(i), apsp, apspmat, apspmap, apspmapback, dstravel, stpaths, dmax, clusterap, changedcluster);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				
			}
			
			//printSuperTargets(sts);
			
			
			return sts;
			
				
			
	}
	
	
	/**
	 * dominated targets are considered as single super target/ cluster
	 * @param tmpgraphmaps
	 * @param attackhistory
	 * @param apsp
	 * @param apspmat
	 * @param apspmap
	 * @param apspmapback
	 * @param dstravel
	 * @param stpaths
	 * @param dmax
	 * @param attackclustershisotry
	 * @param clusteredtargets
	 * @param clustercenters
	 * @param newtargetstocluster
	 * @param domindatednodes
	 * @param RADIUS
	 * @param tmpgraph
	 * @param clusterap
	 * @return
	 */
	private static HashMap<Integer, SuperTarget> constructSuperTargetsV3(HashMap<Integer, TargetNode> tmpgraphmaps,
			 ArrayList<Integer> attackhistory, AllPairShortestPath apsp, int[][] apspmat,
			HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback,
			HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax,
			HashMap<Integer, ArrayList<Integer>> attackclustershisotry, ArrayList<Integer> clusteredtargets, 
			HashMap<Integer,Integer> clustercenters, ArrayList<Integer> newtargetstocluster, ArrayList<TargetNode> domindatednodes, 
			int RADIUS, ArrayList<TargetNode> tmpgraph, HashMap<Integer, int[]> clusterap) {
		
			
			// for now build a cluster
			
			//for(Integer attackt: attackhistory)
			
			//ArrayList<Integer> attackcluster = new ArrayList<Integer>();//createAttackCluster(tmpgraphmaps, attackedtarget);
			
			// base target will alws be alone, so no clustering using base target
		
		
		
		

		/*System.out.println("Current attack clusters ");


		printClusters(attackclustershisotry);

		System.out.println("Current attacked targets ");

		printTargets(newtargetstocluster);


		System.out.println("Clustered targets ");

		printTargets(clusteredtargets);
*/




		System.out.println("Updating attack clusters ");
		
		
		

		ArrayList<Integer> changedcluster = updateClusters(attackclustershisotry, clusteredtargets, clustercenters, newtargetstocluster, tmpgraphmaps, RADIUS, domindatednodes, tmpgraph);







		// create clusters

		int totalcluster = tmpgraphmaps.size() - clusteredtargets.size() + attackclustershisotry.size() /*- domindatednodes.size()*/;


		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[totalcluster];

		for(int i=0; i<totalcluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}

		clusters[0].add(0);

		int clusindex = 1;
		for(ArrayList<Integer> attackcluster: attackclustershisotry.values())
		{

			for(Integer n: attackcluster)
			{

				clusters[clusindex].add(n);
			}
			clusindex++;
		}

		//int index = 2;
		for(Integer t: tmpgraphmaps.keySet())
		{
			TargetNode tnode = tmpgraphmaps.get(t);
			
			if(t != 0 && !clusteredtargets.contains(t) /*&& (!domindatednodes.contains(tnode))*/)
			{
				clusters[clusindex++].add(t);

			}
		}

			
			
			
			//printClusters(clusters);
			
			
			
			

			
			HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,tmpgraphmaps);
			
			//printSuperTargets(sts);
			
			
			System.out.println("\n \n \n After choosing AP \n \n \n");
			
			
			for(int i = 0; i<sts.size(); i++)
			{
				SuperTarget st = sts.get(i);
				
				
				
				
				// if cluster is chnaged then compute AP
				// if not changed then use the previous AP
				
				try {
					chooseDegreeBasedAPV2(sts, tmpgraphmaps, sts.get(i), apsp, apspmat, apspmap, apspmapback, dstravel, stpaths, dmax, clusterap, changedcluster);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				
			}
			
			//printSuperTargets(sts);
			
			
			return sts;
			
				
			
	}
	
	
	/**
	 * it splits the cluster which could not be added to the path
	 * @param tmpgraphmaps
	 * @param attackhistory
	 * @param apsp
	 * @param apspmat
	 * @param apspmap
	 * @param apspmapback
	 * @param dstravel
	 * @param stpaths
	 * @param dmax
	 * @param attackclustershisotry
	 * @param clusteredtargets
	 * @param clustercenters
	 * @param newtargetstocluster
	 * @param domindatednodes
	 * @param RADIUS
	 * @param tmpgraph
	 * @param clusterap
	 * @param notaddedst 
	 * @return
	 */
	private static HashMap<Integer, SuperTarget> constructSuperTargetsV2(HashMap<Integer, TargetNode> tmpgraphmaps,
			 ArrayList<Integer> attackhistory, AllPairShortestPath apsp, int[][] apspmat,
			HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback,
			HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax,
			HashMap<Integer, ArrayList<Integer>> attackclustershisotry, ArrayList<Integer> clusteredtargets, 
			HashMap<Integer,Integer> clustercenters, ArrayList<Integer> newtargetstocluster, ArrayList<TargetNode> domindatednodes, 
			int RADIUS, ArrayList<TargetNode> tmpgraph, HashMap<Integer, int[]> clusterap, ArrayList<Integer> notaddedst) {
		
			
			// for now build a cluster
			
			//for(Integer attackt: attackhistory)
			
			//ArrayList<Integer> attackcluster = new ArrayList<Integer>();//createAttackCluster(tmpgraphmaps, attackedtarget);
			
			// base target will alws be alone, so no clustering using base target
		
		
		
		

		/*System.out.println("Current attack clusters ");


		printClusters(attackclustershisotry);

		System.out.println("Current attacked targets ");

		printTargets(newtargetstocluster);


		System.out.println("Clustered targets ");

		printTargets(clusteredtargets);
*/




		System.out.println("Updating attack clusters ");
		
		
		
		// change history if notaddedst is not empty
		/**
		 * for each ST :
		 * 
		 * 1. split the ST into 2
		 * 2. Update the attack history
		 * 3. update the center
		 */
		
		
		/**
		 * 
		 */
		
		for(Integer sttosplit: notaddedst)
		{
			
			
			sttosplit = sttosplit -1;
			ArrayList<Integer> st = attackclustershisotry.get(sttosplit);

			if(st.size()>1)
			{
				// remove the st
				attackclustershisotry.remove(sttosplit); // give the index/key to remove an object
				int centerid = clustercenters.get(sttosplit);
				//clustercenters.remove(sttosplit);


				// for one cluster, remove the furthest targets from center
				// how many?

				int ntargetstoremove = st.size()/2;

				System.out.println("Targets to remove "+ ntargetstoremove);

				TargetNode src = tmpgraphmaps.get(centerid);
				double[][] ranks = new double[st.size()-1][2];

				int index = 0;

				for(Integer t: st)
				{
					if(t != centerid)
					{
						TargetNode dest = tmpgraphmaps.get(t);
						
						if(src.getNeighbors().contains(dest))
						{
						
						ranks[index][0] = t;
						ranks[index][1] = src.getDistance(dest);
						index++;
						}
						else
						{
							ranks[index][0] = t;
							ranks[index][1] = 9999;
							index++;
						}
					}

				}
				
				
				double[] swap = {0,0};

				for (int k = 0; k < ranks.length; k++) 
				{
					for (int d = 1; d < ranks.length-k; d++) 
					{
						if (ranks[d-1][1] < ranks[d][1])    //descending
						{
							swap = ranks[d];
							ranks[d]  = ranks[d-1];
							ranks[d-1] = swap;
						}
					}
				}
				
				System.out.println("Sorted "+ ntargetstoremove);
				
				ArrayList<Integer> removedtargets = new ArrayList<Integer>();
				
				for(int i=0; i<ntargetstoremove; i++)
				{
					st.remove(st.indexOf((int)ranks[i][0]));
					removedtargets.add((int)ranks[i][0]);
				}
				
				System.out.println("removed "+ ntargetstoremove);
				
				
				
				attackclustershisotry.put(sttosplit, st);
				
				
				
				int newcenterid = removedtargets.get(0);
				
				int newclusterid = attackclustershisotry.size();
				
				attackclustershisotry.put(newclusterid, removedtargets);
				clustercenters.put(newclusterid, newcenterid);
				
				System.out.println("done "+ ntargetstoremove);
				


			}
			
			
			
			
			
			
		}
		
		
		

		ArrayList<Integer> changedcluster = updateClusters(attackclustershisotry, clusteredtargets, clustercenters, newtargetstocluster, tmpgraphmaps, RADIUS, domindatednodes, tmpgraph);







		// create clusters

		int totalcluster = tmpgraphmaps.size() - clusteredtargets.size() + attackclustershisotry.size() - domindatednodes.size();


		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[totalcluster];

		for(int i=0; i<totalcluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}

		clusters[0].add(0);

		int clusindex = 1;
		for(ArrayList<Integer> attackcluster: attackclustershisotry.values())
		{

			for(Integer n: attackcluster)
			{

				clusters[clusindex].add(n);
			}
			clusindex++;
		}

		//int index = 2;
		for(Integer t: tmpgraphmaps.keySet())
		{
			TargetNode tnode = tmpgraphmaps.get(t);
			
			if(t != 0 && !clusteredtargets.contains(t) && (!domindatednodes.contains(tnode)))
			{
				clusters[clusindex++].add(t);

			}
		}

			
			
			
			printClusters(clusters);
			
			
			
			

			
			HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,tmpgraphmaps);
			
			printSuperTargets(sts);
			
			
			System.out.println("\n \n \n After choosing AP \n \n \n");
			
			
			for(int i = 0; i<sts.size(); i++)
			{
				SuperTarget st = sts.get(i);
				
				
				
				
				// if cluster is chnaged then compute AP
				// if not changed then use the previous AP
				
				try {
					chooseDegreeBasedAPV2(sts, tmpgraphmaps, sts.get(i), apsp, apspmat, apspmap, apspmapback, dstravel, stpaths, dmax, clusterap, changedcluster);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				
			}
			
		//	printSuperTargets(sts);
			
			
			return sts;
			
				
			
	}


	private static void printTargets(ArrayList<Integer> targets) {
		
		
		for(Integer t: targets)
		{
			System.out.print(t+ ", ");
		}
		
		System.out.println();
		
	}


	private static ArrayList<Integer> updateClusters(HashMap<Integer, ArrayList<Integer>> attackclustershisotry,
			ArrayList<Integer> clusteredtargets, HashMap<Integer, Integer> clustercenters,
			ArrayList<Integer> newtargetstocluster, HashMap<Integer,TargetNode> tmpgraphmaps, 
			int RADIUS, ArrayList<TargetNode> domindatednodes, ArrayList<TargetNode> tmpgraph) {
		
		
		boolean clusterhappened = false;
		
		ArrayList<Integer> changedcluster = new ArrayList<Integer>();
		
		
		//printNodesWithNeighborsAndPath(tmpgraphmaps);
		
		
		
		if(needToUpdateCluster(newtargetstocluster, clusteredtargets))
		{


			if(newtargetstocluster.contains(0))
			{
				int index = newtargetstocluster.indexOf(0);
				newtargetstocluster.remove(index);
			}


			// check if the cluster history is empty
			if(attackclustershisotry.size()==0)
			{

				for(int tid: newtargetstocluster)
				{
					TargetNode t = tmpgraphmaps.get(tid);


					if(tmpgraph.contains(t))
					{
						//create one cluster with one target from current attacked targets, and that will be the center
						ArrayList<Integer> tmpclus = new ArrayList<Integer>();

						//int tid = currentattackedtargets.get(0);

						tmpclus.add(tid);
						int clusid = attackclustershisotry.size();

						// try to add more targets in the cluster which are not already clustered and not base

						//addMoreToCluster(tmpclus,tid, clusteredtargets, tmpgraphmaps);


						attackclustershisotry.put(clusid, tmpclus);
						clustercenters.put(clusid, tid);
						clusteredtargets.add(tid);
						clusterhappened = true;
						changedcluster.add(clusid);
						break;
					}
				}
			}



			// for every current attacked target
			for(Integer at: newtargetstocluster)
			{
				// if its not clustered already
				TargetNode t = tmpgraphmaps.get(at);
				
				if(!clusteredtargets.contains(at) && tmpgraph.contains(t))
				{
					clusterhappened = true;
					// try to insert it in a cluster with min dist from center
					int assignedcluster = assignToCluster(at, attackclustershisotry, clustercenters, tmpgraphmaps, RADIUS);
					// if -1 then create a new cluster ...
					if(assignedcluster != -1)
					{
						// append at to the assigned cluster
						System.out.println("attckedtarget "+ at + " is assigned to cluster "+ assignedcluster);
						ArrayList<Integer> tmpclus = attackclustershisotry.get(assignedcluster);
						//attackclustershisotry.remove(assignedcluster);
						tmpclus.add(at);
						attackclustershisotry.put(assignedcluster, tmpclus);
						
						if(!changedcluster.contains(assignedcluster))// if it already does not contain the chnaged cluster
						{
							changedcluster.add(assignedcluster);
						}

						// try to add more targets


						//update clustered targets
						clusteredtargets.add(at);

					}
					else if(assignedcluster == -1)
					{
						// create new cluster
						// update center
						// update clustered targets



						ArrayList<Integer> tmpclus = new ArrayList<Integer>();
						tmpclus.add(at);
						int clusid = attackclustershisotry.size();
						attackclustershisotry.put(clusid, tmpclus);
						clustercenters.put(clusid, at);
						clusteredtargets.add(at);
						changedcluster.add(clusid);

						System.out.println("attckedtarget "+ at + " is creating a new cluster with id "+ clusid);


					}
				}
			}

		}

		
		System.out.println("trying to add more targets in the clusters");
		
		int nclus = attackclustershisotry.size();
		
		Iterator<Integer> itr = attackclustershisotry.keySet().iterator();
		
		
		while(itr.hasNext())
		{
			
			int clusid = itr.next();
			
			ArrayList<Integer> tmpclus = attackclustershisotry.get(clusid);
			
			int prevclustersize = tmpclus.size();
			
			int centerid = clustercenters.get(clusid);
			
			//if(tmpclus.size()<5)
			{

				addMoreToCluster(tmpclus, centerid, clusteredtargets, tmpgraphmaps, clusid,RADIUS, tmpgraph);
			}
			
			//attackclustershisotry.remove(clusid);
			
			int currentsize = tmpclus.size();
			
			attackclustershisotry.put(clusid, tmpclus);
			
			if(currentsize>prevclustersize && (!changedcluster.contains(clusid))) // if cluster size increases then add the cluster to changed list
			{
				changedcluster.add(clusid);
			}
			
			
			
		}
		
		
		return changedcluster;
		
	}


	private static boolean needToUpdateCluster(ArrayList<Integer> currentattackedtargets,
			ArrayList<Integer> clusteredtargets) {

		
		
		
		for(Integer t: currentattackedtargets)
		{
			if(!clusteredtargets.contains(t))
				return true;
		}
		
		
		return false;
	}


	private static void addMoreToCluster(ArrayList<Integer> tmpclus, int targetid, ArrayList<Integer> clusteredtargets,
			HashMap<Integer, TargetNode> tmpgraphmaps, int clusid, int RADIUS, ArrayList<TargetNode> tmpgraph) {
		
		
		TargetNode centernode = tmpgraphmaps.get(targetid);
		
		for(TargetNode nei: centernode.getNeighbors())
		{
			if((tmpclus.size()<5) && (nei.getTargetid() != 0) && (!clusteredtargets.contains(nei.getTargetid()))  && (centernode.getDistance(nei) <= RADIUS) && tmpgraph.contains(nei) )
			{
				
				System.out.println("adding target "+ nei.getTargetid() + " to cluster "+ clusid + " dist("+centernode.getTargetid()+","+nei.getTargetid()+ "="+centernode.getDistance(nei));
				tmpclus.add(nei.getTargetid());
				clusteredtargets.add(nei.getTargetid());
			}
		}
		
	}


	private static int assignToCluster(Integer attackedtarget, HashMap<Integer, ArrayList<Integer>> attackclustershisotry,
			HashMap<Integer, Integer> clustercenters, HashMap<Integer, TargetNode> tmpgraphmaps, int RADIUS) {
		
		
		double mindist = Double.MAX_VALUE;
		int minclus = -1;
		
		
		for(Integer clusid : attackclustershisotry.keySet())
		{
			ArrayList<Integer> clus = attackclustershisotry.get(clusid);
			int center = clustercenters.get(clusid);
			
			if(clus.size()==5)
				continue;
			
			// find the distance from center to attackedtarget
			TargetNode centertarget = tmpgraphmaps.get(center);
			TargetNode attackednode = tmpgraphmaps.get(attackedtarget);
			
			// if they are neighbor then compute distance
			
			if(centertarget.getNeighbors().contains(attackednode))
			{
				double dist = centertarget.getDistance(attackednode);
				if((dist<=RADIUS) && (dist<mindist) && (clus.size()<5))
				{
					mindist = dist;
					minclus = clusid;
				}
			}
			// else move to next cluster
			
		}
		
		return minclus;
	}


	private static void chooseAP(HashMap<Integer, SuperTarget> sts, HashMap<Integer, TargetNode> targetmaps, SuperTarget tempst,
			 AllPairShortestPath apsp,
				int[][] apspmat, HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback, 
				HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax, HashMap<Integer,int[]> clusterap, ArrayList<Integer> changedcluster) {
		
		//for(SuperTarget tempst: sts.values())
		{
			//SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
			//tempst.stid = 200 + st1.stid + st2.stid;
			// printSuperTarget(tempst);
			/**
			 * for every pair of access points
			 */
			
			
			
			
			double mindi = Double.MAX_VALUE;
			int stid1 = -1;
			//int stid2 = -1;
			int aid1 = -1;
			int aid2 = -1;
			int s1id = -1;
			int s2id = -2;
			int s1ap = -1;
			int s2ap = -1;
			double sda1a2 = -1;
			double sd1 = -1;
			double sd2 = -1;
			ArrayList<Integer> a1a2path = new ArrayList<Integer>();
			ArrayList<Integer> a2path = new ArrayList<Integer>();
			ArrayList<Integer> a1path = new ArrayList<Integer>();
			SuperTarget dminst = new SuperTarget();
			
			double mindist = 500000;
			
			ArrayList<int[]> doneappair = new ArrayList<int[]>();
			ArrayList<int[]> donestpair = new ArrayList<int[]>();
			
			
			
			

			if(tempst.nodes.size()>1)
			{
				
				System.out.println("Computing ap for ST "+ tempst.stid);
				
				if(tempst.stid==15)
				{
					System.out.println("shortestdist(a1,s1) ");
				}



				//check if 0 is its neighbor if so just keep the ap that connects to base

				if(tempst.neighbors.containsKey(0) && sts.size()>2)
				{
					//compute every parameter

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								// try compute travelling path for cluster
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								// for every pair of supertargets which are not st1 or st2
								
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if((s1.stid != s2.stid) && (notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else if(tempst.neighbors.containsKey(0) && sts.size()==2)
				{
					//compute every parameter

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								
								// for every pair of supertargets which are not st1 or st2
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if(/*(s1.stid != s2.stid) && */(notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);
									


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											//ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2


											//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
											//double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
											

											//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);



											//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

											/*if(dista1a2 == 0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

												dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
											}*/

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else
				{
					
					System.out.println("Here for ST" + tempst.stid);
					
					
					//printSuperTarget(tempst);

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							boolean notd = notDone(a1.getTargetid(), a2.getTargetid(), doneappair);
							if((a1.getTargetid() != a2.getTargetid()) && (notd))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								// for every pair of supertargets which are not st1 or st2
								for(SuperTarget s1: tempst.neighbors.values()) // need neighbor cluster of a1
								{
									for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
									{

										if(/*(s1.stid != s2.stid) && */(notDone(s1.stid, s2.stid, donestpair)))
										{
											
											int[] stpair = {s1.stid, s2.stid};
											donestpair.add(stpair);
										


											boolean arebothnei = areBothNei(s1,s2,tempst);

											if(s1.stid != tempst.stid && 
													s2.stid != tempst.stid  && arebothnei)
											{
												// measure the distance between a1->s1 and a2->s2

												//a1->s1 distance between a target and a supertarget
												/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


												//ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
												ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
												ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


												double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
												//TODO
												//System.out.println("shortestdist(a1,s1) "+ d1);

												double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
												//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
												//TODO
												//System.out.println("shortestdist(a2,s2) "+ d2);

												// next measure the intra cluster shortest traveling path using a1 and a2


												//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
												
												/*double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
												
												
												
												//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);



												//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

												if(dista1a2 == 0)
													//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												{
													//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

													dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
												}
*/
												//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

												// if any of the dist is <0 we know that it's not possible to have a path

												if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
												{
													double totaldi = d1[0] + dista1a2 + d2[0];
													if(totaldi < mindi)
													{
														//stid1 = .stid;
														//stid2 = st2.stid;
														dminst = tempst;
														aid1 = a1.getTargetid();
														aid2 = a2.getTargetid();
														mindi = totaldi;
														sda1a2 = dista1a2;
														s1id = s1.stid;
														s2id = s2.stid;
														sd1 = d1[0];
														sd2= d2[0];
														s1ap = (int)d1[1];
														s2ap = (int)d2[1];
														//spath.add(tmpspath.get(0));
														a1a2path.clear();
														for(Integer in: tmpa1a2spath)
														{
															a1a2path.add(in);
														}


														a1path.clear();
														for(Integer in: tmpa1path)
														{
															a1path.add(in);
														}


														a2path.clear();
														for(Integer in: tmpa2path)
														{
															a2path.add(in);
														}

														/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


													}
												}

											}
										}
									}
								}
							}
						}
					}
				}
			}

			//}
			//						}
			//
			//					}
			//				}
			//			}


			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
			//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			/**
			 * Merge two cluster which have min d
			 * For the new id of supertargets concat the strings with comma
			 * remove the old two clusters
			 * add the new one. 
			 */
			
			
			//TODO what to do when there is no path from a1 to a2?????
			
			
			
			if(sda1a2==-1 && tempst.nodes.size()>1)
			{
				System.out.println("shortestdist(a1,s1) ");
				//printSuperTarget(tempst);
			}

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				
				// Issue: How to fix the issue : connection to base node might be cut
				// Fix: Always keep the ap which connects to base node
				
				
				
				System.out.println("AP done for st "+ tempst.stid);
				
				dstravel.put(tempst.stid, sda1a2);
				stpaths.put(tempst.stid, a1a2path);
				//SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
				// printSuperTarget(newst);
				/*sts.remove(stid1);
				sts.remove(stid2);*/
				//printSuperTargets(sts);
				
				updateAP(tempst, sts, aid1, aid2);
				
				
				clusterap.put(tempst.stid, new int[] {aid1, aid2});
				
				//sts.put(newst.stid, newst);
				//update the neighbors of ST
				// System.out.println("\n\n After merging # supertargets : "+ sts.size());
				// System.out.println("\n After merging new supertargets : ");
				
				/**
				 * add s1 and s2 as neighbor of st if they are not already
				 */
				
				addNei(tempst.stid, aid1, s1id, sd1, a1path, sts, targetmaps, s1ap);
				addNei(tempst.stid, aid2, s2id, sd2, a2path, sts, targetmaps, s2ap);
				//System.out.println("stid "+ tempst.stid);
				addNei(aid1,aid2,sda1a2, a1a2path, tempst, targetmaps);
				
				
				/**
				 * add a1 a2 as neighbor if they are not already
				 */
				
				
				//printSuperTargets(sts);
				//printNodesWithNeighborsAndPath(targetmaps);
				
				
				System.out.println("a1 "+ aid1 + " a2 "+ aid2 + " for cluster  "+ tempst.stid);
				 
			}
			
		}

		
		
	}
	
	
	
	private static void chooseDegreeBasedAP(HashMap<Integer, SuperTarget> sts, HashMap<Integer, TargetNode> targetmaps, SuperTarget tempst,
			 AllPairShortestPath apsp,
				int[][] apspmat, HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback, 
				HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax) throws Exception {
		
		//for(SuperTarget tempst: sts.values())
		{
			//SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
			//tempst.stid = 200 + st1.stid + st2.stid;
			// printSuperTarget(tempst);
			/**
			 * for every pair of access points
			 */
			
			
			
			
			double mindi = Double.MAX_VALUE;
			int stid1 = -1;
			//int stid2 = -1;
			int aid1 = -1;
			int aid2 = -1;
			int s1id = -1;
			int s2id = -2;
			int s1ap = -1;
			int s2ap = -1;
			double sda1a2 = -1;
			double sd1 = -1;
			double sd2 = -1;
			ArrayList<Integer> a1a2path = new ArrayList<Integer>();
			ArrayList<Integer> a2path = new ArrayList<Integer>();
			ArrayList<Integer> a1path = new ArrayList<Integer>();
			SuperTarget dminst = new SuperTarget();
			
			double mindist = 500000;
			
			ArrayList<int[]> doneappair = new ArrayList<int[]>();
			ArrayList<int[]> donestpair = new ArrayList<int[]>();
			
			int category = -1;
			
			
			
			

			if(tempst.nodes.size()>1)
			{
				
				System.out.println("Computing ap for ST "+ tempst.stid);
				
				if(tempst.stid==15)
				{
					System.out.println("shortestdist(a1,s1) ");
				}



				//check if 0 is its neighbor if so just keep the ap that connects to base

				if(tempst.neighbors.containsKey(0) && sts.size()>2)
				{
					//compute every parameter
					
					category = 1;
					System.out.println("category "+ category);

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								// try compute travelling path for cluster
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0 || dista1a2>dmax)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								// for every pair of supertargets which are not st1 or st2
								
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if((s1.stid != s2.stid) && (notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else if(tempst.neighbors.containsKey(0) && sts.size()==2)
				{
					//compute every parameter
					
					
					category = 2;
					System.out.println("category "+ category);

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0 || dista1a2>dmax)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								
								// for every pair of supertargets which are not st1 or st2
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if(/*(s1.stid != s2.stid) && */(notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);
									


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											//ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2


											//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
											//double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
											

											//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);



											//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

											/*if(dista1a2 == 0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

												dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
											}*/

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else
				{

					System.out.println("Here for ST" + tempst.stid);

					category = 3;
					System.out.println("category "+ category);

					// rank the APs based on their degree

					int [][] rankedap = rankAPBasedOnDegree(tempst, sts);

					int a1 = rankedap[0][0];
					int a2 = rankedap[1][0];




					// next measure the intra cluster shortest traveling path using a1 and a2
					ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
					//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1, a2, tmpa1a2spath);
					//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
					//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					if(dista1a2 == 0 || dista1a2 > dmax)
					{
						//throw new Exception("No path found to compute AP for st "+ tempst.stid);
						//System.out.println("Disjpint cluster need longer path "+ tempst.stid);
						TargetNode a1node = tempst.nodes.get(a1);
						TargetNode a2node = tempst.nodes.get(a2);
						dista1a2 = shortestdist(a1node,a2node, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
					}



					if(/*d1[0] > 0 && d2[0] > 0 &&*/ dista1a2 > 0)
					{
						double totaldi =  dista1a2;
						if(totaldi < mindi)
						{
							//stid1 = .stid;
							//stid2 = st2.stid;
							dminst = tempst;
							aid1 = a1;
							aid2 = a2;
							mindi = totaldi;
							sda1a2 = dista1a2;
							a1a2path.clear();
							for(Integer in: tmpa1a2spath)
							{
								a1a2path.add(in);
							}



						}
					}

				}
			}

			
			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
			//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			/**
			 * Merge two cluster which have min d
			 * For the new id of supertargets concat the strings with comma
			 * remove the old two clusters
			 * add the new one. 
			 */
			
			
			//TODO what to do when there is no path from a1 to a2?????
			
			
			
			if(sda1a2==-1 && tempst.nodes.size()>1)
			{
				System.out.println("shortestdist(a1,s1) ");
				//printSuperTarget(tempst);
			}

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				
				// Issue: How to fix the issue : connection to base node might be cut
				// Fix: Always keep the ap which connects to base node
				
				
					System.out.println("AP done for st "+ tempst.stid);
					
					dstravel.put(tempst.stid, sda1a2);
					stpaths.put(tempst.stid, a1a2path);
					//SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
					// printSuperTarget(newst);
					/*sts.remove(stid1);
					sts.remove(stid2);*/
					//printSuperTargets(sts);
					
					updateAP(tempst, sts, aid1, aid2);
					
					
				//	clusterap.put(tempst.stid, new int[] {aid1, aid2});
					
					//sts.put(newst.stid, newst);
					//update the neighbors of ST
					// System.out.println("\n\n After merging # supertargets : "+ sts.size());
					// System.out.println("\n After merging new supertargets : ");
					
					/**
					 * add s1 and s2 as neighbor of st if they are not already
					 */
					
					if(category != 3)
					{
					
						addNei(tempst.stid, aid1, s1id, sd1, a1path, sts, targetmaps, s1ap);
						addNei(tempst.stid, aid2, s2id, sd2, a2path, sts, targetmaps, s2ap);
					}
					//System.out.println("stid "+ tempst.stid);
					addNei(aid1,aid2,sda1a2, a1a2path, tempst, targetmaps);
					
					
					/**
					 * add a1 a2 as neighbor if they are not already
					 */
					
					
					//printSuperTargets(sts);
					//printNodesWithNeighborsAndPath(targetmaps);
					
					
					System.out.println("a1 "+ aid1 + " a2 "+ aid2 + " for cluster  "+ tempst.stid);
				
				
				
					 
			}
			
		}

		
		
	}
	
	
	private static void chooseDegreeBasedAP(HashMap<Integer, SuperTarget> sts, HashMap<Integer, TargetNode> targetmaps, SuperTarget tempst,
			 AllPairShortestPath apsp,
				int[][] apspmat, HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback, 
				HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax, HashMap<Integer,int[]> clusterap, ArrayList<Integer> changedcluster) throws Exception {
		
		//for(SuperTarget tempst: sts.values())
		{
			//SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
			//tempst.stid = 200 + st1.stid + st2.stid;
			// printSuperTarget(tempst);
			/**
			 * for every pair of access points
			 */
			
			
			
			
			double mindi = Double.MAX_VALUE;
			int stid1 = -1;
			//int stid2 = -1;
			int aid1 = -1;
			int aid2 = -1;
			int s1id = -1;
			int s2id = -2;
			int s1ap = -1;
			int s2ap = -1;
			double sda1a2 = -1;
			double sd1 = -1;
			double sd2 = -1;
			ArrayList<Integer> a1a2path = new ArrayList<Integer>();
			ArrayList<Integer> a2path = new ArrayList<Integer>();
			ArrayList<Integer> a1path = new ArrayList<Integer>();
			SuperTarget dminst = new SuperTarget();
			
			double mindist = 500000;
			
			ArrayList<int[]> doneappair = new ArrayList<int[]>();
			ArrayList<int[]> donestpair = new ArrayList<int[]>();
			
			int category = -1;
			
			
			
			

			if(tempst.nodes.size()>1)
			{
				
				System.out.println("Computing ap for ST "+ tempst.stid);
				
				if(tempst.stid==5)
				{
					System.out.println("shortestdist(a1,s1) ");
				}



				//check if 0 is its neighbor if so just keep the ap that connects to base

				if(tempst.neighbors.containsKey(0) && sts.size()>2)
				{
					//compute every parameter
					
					category = 1;
					System.out.println("category "+ category);

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								// try compute travelling path for cluster
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								// for every pair of supertargets which are not st1 or st2
								
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if((s1.stid != s2.stid) && (notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else if(tempst.neighbors.containsKey(0) && sts.size()==2)
				{
					//compute every parameter
					
					
					category = 2;
					System.out.println("category "+ category);

					/*int aid1 = -1;
					int aid2 = -1;
					int s1id = -1;
					int s2id = -2;
					int s1ap = -1;
					int s2ap = -1;
					double sda1a2 = -1;
					double sd1 = -1;
					double sd2 = -1;
					ArrayList<Integer> a1a2path = new ArrayList<Integer>();
					ArrayList<Integer> a2path = new ArrayList<Integer>();
					ArrayList<Integer> a1path = new ArrayList<Integer>();
					 */

					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								
								// for every pair of supertargets which are not st1 or st2
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if(/*(s1.stid != s2.stid) && */(notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);
									


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											//ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2


											//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
											//double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
											

											//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);



											//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

											/*if(dista1a2 == 0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

												dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
											}*/

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else
				{

					System.out.println("Here for ST" + tempst.stid);

					category = 3;
					System.out.println("category "+ category);

					// rank the APs based on their degree

					int [][] rankedap = rankAPBasedOnDegree(tempst,sts);

					int a1 = rankedap[0][0];
					int a2 = rankedap[1][0];




					// next measure the intra cluster shortest traveling path using a1 and a2
					ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
					//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1, a2, tmpa1a2spath);
					//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
					//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					if(dista1a2 == 0)
					{
						throw new Exception("No path found to compute AP for st "+ tempst.stid);
						//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

						//dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
					}



					if(/*d1[0] > 0 && d2[0] > 0 &&*/ dista1a2 > 0)
					{
						double totaldi =  dista1a2;
						if(totaldi < mindi)
						{
							//stid1 = .stid;
							//stid2 = st2.stid;
							dminst = tempst;
							aid1 = a1;
							aid2 = a2;
							mindi = totaldi;
							sda1a2 = dista1a2;
							a1a2path.clear();
							for(Integer in: tmpa1a2spath)
							{
								a1a2path.add(in);
							}



						}
					}

				}
			}

			
			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
			//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			/**
			 * Merge two cluster which have min d
			 * For the new id of supertargets concat the strings with comma
			 * remove the old two clusters
			 * add the new one. 
			 */
			
			
			//TODO what to do when there is no path from a1 to a2?????
			
			
			
			if(sda1a2==-1 && tempst.nodes.size()>1)
			{
				System.out.println("shortestdist(a1,s1) ");
				//printSuperTarget(tempst);
			}

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				
				// Issue: How to fix the issue : connection to base node might be cut
				// Fix: Always keep the ap which connects to base node
				
				
					System.out.println("AP done for st "+ tempst.stid);
					
					dstravel.put(tempst.stid, sda1a2);
					stpaths.put(tempst.stid, a1a2path);
					//SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
					// printSuperTarget(newst);
					/*sts.remove(stid1);
					sts.remove(stid2);*/
					//printSuperTargets(sts);
					
					updateAP(tempst, sts, aid1, aid2);
					
					
					clusterap.put(tempst.stid, new int[] {aid1, aid2});
					
					//sts.put(newst.stid, newst);
					//update the neighbors of ST
					// System.out.println("\n\n After merging # supertargets : "+ sts.size());
					// System.out.println("\n After merging new supertargets : ");
					
					/**
					 * add s1 and s2 as neighbor of st if they are not already
					 */
					
					if(category != 3)
					{
					
						addNei(tempst.stid, aid1, s1id, sd1, a1path, sts, targetmaps, s1ap);
						addNei(tempst.stid, aid2, s2id, sd2, a2path, sts, targetmaps, s2ap);
					}
					//System.out.println("stid "+ tempst.stid);
					addNei(aid1,aid2,sda1a2, a1a2path, tempst, targetmaps);
					
					
					/**
					 * add a1 a2 as neighbor if they are not already
					 */
					
					
					//printSuperTargets(sts);
					//printNodesWithNeighborsAndPath(targetmaps);
					
					
					System.out.println("a1 "+ aid1 + " a2 "+ aid2 + " for cluster  "+ tempst.stid);
				
				
				
					 
			}
			
		}

		
		
	}
	
	
	/**
	 * choose ap based on degree of node for all of the clusters where cluster size > 2
	 * 
	 * @param sts
	 * @param targetmaps
	 * @param tempst
	 * @param apsp
	 * @param apspmat
	 * @param apspmap
	 * @param apspmapback
	 * @param dstravel
	 * @param stpaths
	 * @param dmax
	 * @param clusterap
	 * @param changedcluster
	 * @throws Exception
	 */
	private static void chooseDegreeBasedAPV2(HashMap<Integer, SuperTarget> sts, HashMap<Integer, TargetNode> targetmaps, SuperTarget tempst,
			 AllPairShortestPath apsp,
				int[][] apspmat, HashMap<Integer,Integer> apspmap, HashMap<Integer,Integer> apspmapback, 
				HashMap<Integer,Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, int dmax, HashMap<Integer,int[]> clusterap, ArrayList<Integer> changedcluster) throws Exception {
		
		//for(SuperTarget tempst: sts.values())
		{
			//SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
			//tempst.stid = 200 + st1.stid + st2.stid;
			// printSuperTarget(tempst);
			/**
			 * for every pair of access points
			 */
			
			
			
			
			double mindi = Double.MAX_VALUE;
			int stid1 = -1;
			//int stid2 = -1;
			int aid1 = -1;
			int aid2 = -1;
			int s1id = -1;
			int s2id = -2;
			int s1ap = -1;
			int s2ap = -1;
			double sda1a2 = -1;
			double sd1 = -1;
			double sd2 = -1;
			ArrayList<Integer> a1a2path = new ArrayList<Integer>();
			ArrayList<Integer> a2path = new ArrayList<Integer>();
			ArrayList<Integer> a1path = new ArrayList<Integer>();
			SuperTarget dminst = new SuperTarget();
			
			double mindist = 500000;
			
			ArrayList<int[]> doneappair = new ArrayList<int[]>();
			ArrayList<int[]> donestpair = new ArrayList<int[]>();
			
			int category = -1;
			
			
			
			

			if(tempst.nodes.size()>1)
			{
				
				System.out.println("Computing ap for ST "+ tempst.stid);
				
				

				//check if 0 is its neighbor if so just keep the ap that connects to base

				if(tempst.neighbors.containsKey(0) && sts.size()>2)
				{
					System.out.println("Here for ST" + tempst.stid);

					category = 1;
					System.out.println("category "+ category);

					// rank the APs based on their degree

					int [][] rankedap = rankAPBasedOnDegree(tempst,sts);

					int a1 = rankedap[0][0];
					int a2 = getAPWithBase(tempst, a1);//get the AP which is connected with base station and different than a1
					if(a2 == -1)
					{ 
						a2 = rankedap[1][0];
					}




					// next measure the intra cluster shortest traveling path using a1 and a2
					ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
					//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1, a2, tmpa1a2spath);
					//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
					//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					if(dista1a2 == 0 || dista1a2>dmax)
					{
						//throw new Exception("No path found to compute AP for st "+ tempst.stid);
						//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

						TargetNode a1node = tempst.nodes.get(a1);
						TargetNode a2node = tempst.nodes.get(a2);
						
						dista1a2 = shortestdist(a1node,a2node, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
					}



					if(/*d1[0] > 0 && d2[0] > 0 &&*/ dista1a2 > 0)
					{
						double totaldi =  dista1a2;
						if(totaldi < mindi)
						{
							//stid1 = .stid;
							//stid2 = st2.stid;
							dminst = tempst;
							aid1 = a1;
							aid2 = a2;
							mindi = totaldi;
							sda1a2 = dista1a2;
							a1a2path.clear();
							for(Integer in: tmpa1a2spath)
							{
								a1a2path.add(in);
							}



						}
					}


				}
				else if(tempst.neighbors.containsKey(0) && sts.size()==2)
				{
					//compute every parameter
					
					// no need to compute anything
					// just compute intra cluster travelling distance
					
					
					category = 2;
					System.out.println("category "+ category);

					
					for(TargetNode a1: tempst.ap.values())
					{
						/**
						 * 	for every pair of other supertargets which are not st1 and st2
						 * */
						for(TargetNode a2: tempst.ap.values())
						{
							/**
							 * 		find the min d = d1 + dist(a1,a2) + d2
							 */
							// should the ap be same for entry and exit ?
							if((a1.getTargetid() != a2.getTargetid()) && (notDone(a1.getTargetid(), a2.getTargetid(), doneappair)))
							{
								
								if(a1.getTargetid()==224 && a2.getTargetid()==272)
								{
									System.out.println("x");
								}
								
								
								int[] appair = {a1.getTargetid(), a2.getTargetid()};
								doneappair.add(appair);
								
								// next measure the intra cluster shortest traveling path using a1 and a2
								ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
								//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
								//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
								//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
								if(dista1a2 == 0 || dista1a2>dmax)
								{
									//throw new Exception("No path found to compute AP for st "+ tempst.stid);
									//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

									dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
								}
								
								
								
								// for every pair of supertargets which are not st1 or st2
								SuperTarget s1 = sts.get(0); // get base, s1 always will be base

								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{

									if(/*(s1.stid != s2.stid) && */(notDone(s1.stid, s2.stid, donestpair)))
									{
										
										int[] stpair = {s1.stid, s2.stid};
										donestpair.add(stpair);
									


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											//ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();


											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2


											//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
											//double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), tmpa1a2spath);
											

											//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);



											//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

											/*if(dista1a2 == 0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												//System.out.println("Disjpint cluster need longer path "+ tempst.stid);

												dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, (int)dmax); 
											}*/

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}


													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}


													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


												}
											}

										}
									}
								}

							}
						}
					}



				}
				else
				{

					System.out.println("Here for ST" + tempst.stid);

					category = 3;
					System.out.println("category "+ category);

					// rank the APs based on their degree

					int [][] rankedap = rankAPBasedOnDegree(tempst, sts);

					int a1 = rankedap[0][0];
					int a2 = rankedap[1][0];




					// next measure the intra cluster shortest traveling path using a1 and a2
					ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
					//double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					double dista1a2 = travelingGreedyPathAP(tempst.nodes, dmax, tempst.nodes.size(), a1, a2, tmpa1a2spath);
					//double dddd = pathAPForST(tempst.nodes, dmax, tempst.nodes.size(), a1.getTargetid(), a2.getTargetid(), 1, tmpa1a2spath);
					//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
					if(dista1a2 == 0 || dista1a2>dmax)
					{
						//throw new Exception("No path found to compute AP for st "+ tempst.stid);
						//System.out.println("Disjpint cluster need longer path "+ tempst.stid);
						TargetNode a1node = tempst.nodes.get(a1);
						TargetNode a2node = tempst.nodes.get(a2);
						dista1a2 = shortestdist(a1node,a2node, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst, dmax); 
					}



					if(/*d1[0] > 0 && d2[0] > 0 &&*/ dista1a2 > 0)
					{
						double totaldi =  dista1a2;
						if(totaldi < mindi)
						{
							//stid1 = .stid;
							//stid2 = st2.stid;
							dminst = tempst;
							aid1 = a1;
							aid2 = a2;
							mindi = totaldi;
							sda1a2 = dista1a2;
							a1a2path.clear();
							for(Integer in: tmpa1a2spath)
							{
								a1a2path.add(in);
							}



						}
					}

				}
			}

			
			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
			//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			/**
			 * Merge two cluster which have min d
			 * For the new id of supertargets concat the strings with comma
			 * remove the old two clusters
			 * add the new one. 
			 */
			
			
			//TODO what to do when there is no path from a1 to a2?????
			
			
			
			if(sda1a2==-1 && tempst.nodes.size()>1)
			{
				System.out.println("shortestdist(a1,s1) ");
				//printSuperTarget(tempst);
			}

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				
				// Issue: How to fix the issue : connection to base node might be cut
				// Fix: Always keep the ap which connects to base node
				
				
					System.out.println("AP done for st "+ tempst.stid);
					
					dstravel.put(tempst.stid, sda1a2);
					stpaths.put(tempst.stid, a1a2path);
					//SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
					// printSuperTarget(newst);
					/*sts.remove(stid1);
					sts.remove(stid2);*/
					//printSuperTargets(sts);
					
					updateAP(tempst, sts, aid1, aid2);
					
					
					//clusterap.put(tempst.stid, new int[] {aid1, aid2});
					
					//sts.put(newst.stid, newst);
					//update the neighbors of ST
					// System.out.println("\n\n After merging # supertargets : "+ sts.size());
					// System.out.println("\n After merging new supertargets : ");
					
					/**
					 * add s1 and s2 as neighbor of st if they are not already
					 */
					
					if(category == 2)
					{
					
						addNei(tempst.stid, aid1, s1id, sd1, a1path, sts, targetmaps, s1ap);
						addNei(tempst.stid, aid2, s2id, sd2, a2path, sts, targetmaps, s2ap);
					}
					//System.out.println("stid "+ tempst.stid);
					addNei(aid1,aid2,sda1a2, a1a2path, tempst, targetmaps);
					
					
					/**
					 * add a1 a2 as neighbor if they are not already
					 */
					
					
					//printSuperTargets(sts);
					//printNodesWithNeighborsAndPath(targetmaps);
					
					
					System.out.println("a1 "+ aid1 + " a2 "+ aid2 + " for cluster  "+ tempst.stid);
				
				
				
					 
			}
			
		}

		
		
	}
	
	
	
	
	
	private static int getAPWithBase(SuperTarget tempst, int a1) {
		
		
		for(TargetNode t: tempst.ap.values())
		{
			for(TargetNode nei: t.getNeighbors())
			{
				if(nei.getTargetid()==0 && (t.getTargetid() != a1))
				{
					return t.getTargetid();
				}
			}
		}
		
		
		
		return -1;
	}


	private static int [][] rankAPBasedOnDegree(SuperTarget tempst, HashMap<Integer,SuperTarget> sts) {
		
		
		//ArrayList<Integer> rank = new ArrayList<Integer>();
		
		
		int [][] ranks = new int[tempst.ap.size()][2];
		
		int index=0;
		
		for(TargetNode ap: tempst.ap.values())
		{
			int degree = degreeAP(ap, tempst, sts);
			ranks[index][0] = ap.getTargetid();
			ranks[index][1] = degree;
			index++;
		}
		
		int[] swap = {0,0};

		for (int k = 0; k < ranks.length; k++) 
		{
			for (int d = 1; d < ranks.length-k; d++) 
			{
				if (ranks[d-1][1] < ranks[d][1])    // ascending order
				{
					swap = ranks[d];
					ranks[d]  = ranks[d-1];
					ranks[d-1] = swap;
				}
			}
		}
		
		
		return ranks;
	}


	private static int degreeAP(TargetNode ap, SuperTarget tempst, HashMap<Integer,SuperTarget> sts) {
		
		
		
		int count = 0;
		
		for(TargetNode nei: ap.getNeighbors())
		{
			
			for(SuperTarget st: sts.values())
			{
				if((st.stid!=tempst.stid) && st.nodes.keySet().contains(nei.getTargetid()))
				{
					count++;
					
				}
			}
		}
		
		
		return count;
	}


	private static boolean notDone(int targetid, int targetid2, ArrayList<int[]> doneappair) {
		
		
		for(int[] pair: doneappair)
		{

			if((pair[0] == targetid && pair[1] == targetid2) || (pair[0] == targetid2 && pair[1] == targetid) )
			{
				return false;
			}
		}
		
		
		return true;
	}


	public static double travelingGreedyPathAP(HashMap<Integer, TargetNode> targetmaps, double dmax, int nTargets, int base, int dest, ArrayList<Integer> path) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */
		
		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		
		for(TargetNode t: targetmaps.values())
		{
			targets.add(t);
		}

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrixSuperT(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		SecurityGameContraction.purifyAPSPMatrixZero(apspmat, targets, nTargets, map, mapback);
		



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		 
		
		int s = map.get(base);
		int d = map.get(dest);
		

		double dist = travelingPathAP(base, dest, targets, dmax, targetssorted, apspmat, map,mapback, path);
		
		/*for(Integer n: tcur)
		{
			path.add(n);
		}*/

		return dist;
	}
	
	
	public static void makeAdjacencyMatrixSuperT(int[][] adjacencymatrix,
			ArrayList<TargetNode> targets, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				//System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");


				if(targets.contains(nei))
				{
					adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
				}
			}
			i++;
		}


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();
				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}





	}
	
	
	
	public static double travelingPathAP(int base, int dest, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, ArrayList<Integer> path) 
			{
		
		
		
		
		
		
		


		int i=0; 

		ArrayList<Integer> tcur = new ArrayList<Integer>();
		ArrayList<Integer> greedypath = new ArrayList<Integer>();

		if(targetssorted[0][0]==base)
			i=1;
		else 
			i=0;


		int totaldist = 0;
		//for(int res = 0; res<nRes; res++)
		{
			/*System.out.println("res = "+res);
			System.out.println("Adding base to greedy path : ");*/
			greedypath.clear();
			greedypath.add(base);
			greedypath.add(dest);
			//printGreedyPath(greedypath);

			
			if(base!=dest)
			{
				totaldist = apspmat[map.get(base)][map.get(dest)];
			}
			boolean continueFlag= false;
			while ( i<targetssorted.length && continueFlag==false)
			{
				
				int besttotaldisttemp = -1;
				int bestj = -1;
				/*System.out.println("i = "+i);
				System.out.println("continueFlag : "+continueFlag);
				System.out.println("besttotaldisttemp = "+besttotaldisttemp);
				System.out.println("bestj : "+bestj);*/
				
				//System.out.println("trying to insert  : "+targetssorted[i][0]);

				for(int j=1; j< greedypath.size(); j++ )
				{
					int s = map.get(greedypath.get(j-1));
					int d = map.get(greedypath.get(j));
					int totaldisttemp = totaldist - apspmat[s][d];

					s = map.get(greedypath.get(j-1));
					d = map.get(targetssorted[i][0]);
					totaldisttemp +=  apspmat[s][d];


					s = map.get(targetssorted[i][0]);
					d = map.get(greedypath.get(j));
					totaldisttemp +=  apspmat[s][d];

					/*System.out.println("totaldisttemp = "+totaldisttemp);
					System.out.println("bestj : "+bestj);*/

					if (totaldisttemp>0 &&  (totaldisttemp<=dmax) && 
							((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) )
					{
						bestj = j;
						besttotaldisttemp = totaldisttemp;
						/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
						System.out.println("updating bestj : "+bestj);*/
					}
				}
				if(bestj!=-1)
				{
					//System.out.println("could not insert, bestj : "+bestj);
					//continueFlag = true;
					//break;
					//System.out.println("inserting : "+targetssorted[i][0]+" , greedy path : ");
					greedypath = insertTsrt(greedypath, 0, bestj-1, targetssorted[i][0], bestj, greedypath.size()-1);
					//printGreedyPath(greedypath);
					totaldist = besttotaldisttemp;
				}
				
				i=i+1;
			}//while loop
			
			
			
			
			

			/*if(i==targetssorted.length)
			{
				System.out.println("i= "+i+"=targetssorted.length");
				System.out.println("All targets can be added to tcur: ");
				for(int k=0; k<targetssorted.length; k++)
				{
					tcur.add(targetssorted[k][0]);
				}
				//printGreedyPath(tcur);
				return tcur;
			}*/
			//System.out.println("res= "+res+", greedypath : ");
			//printGreedyPath(greedypath);
		}//end of for loop
		
	//	path = greedypath;
		
		for(Integer n: greedypath)
		{
			path.add(n);
		}


		return totaldist;
		
		

		/*int ulimit = targetssorted[i][1];
		System.out.println("ulimit "+ulimit);
		System.out.println("tcur : ");

		
		//if(greedypath.size()==2)
		
			//fidn the min
			// set the ulimit
			
		
			double min = getTargetNode(base, targets).attackerreward;
			for(int p=0; p<greedypath.size(); p++)
			{
			
			if(getTargetNode(greedypath.get(p), targets).attackerreward<min)
				min = (int)getTargetNode(greedypath.get(p), targets).attackerreward;
			
		
			}
		ulimit=(int)min;
		
		for(int k=0; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				tcur.add(targetssorted[k][0]);
			}
		}
		//printGreedyPath(tcur);
		return tcur;*/


			}

	public static ArrayList<Integer> insertTsrt(ArrayList<Integer> greedypath, int i, int j,
			int k, int l, int m) {


		
		
		
		for(int p=0; p<greedypath.size(); p++)
		{
			if(k==greedypath.get(p))
				return greedypath;
		}
		

		ArrayList<Integer> tmp = new ArrayList<Integer>();

		for(int p = i; p<=j; p++)
		{
			tmp.add(greedypath.get(p));
		}

		tmp.add(k);


		for(int p = l; p<=m; p++)
		{
			tmp.add(greedypath.get(p));
		}
		return tmp;


	}
	
	


	private static ArrayList<Integer> createAttackCluster(HashMap<Integer, TargetNode> tmpgraphmaps,
			int attackedtarget) {
		
		
		ArrayList<Integer> attackcluster = new ArrayList<Integer>();
		
		
		if(tmpgraphmaps.keySet().contains(attackedtarget))
		{
			System.out.println("Attack target "+ attackedtarget +" is in tmpraph \n Neighbors are :");
			attackcluster.add(attackedtarget);
			// build cluster around attacked target if they are in tmpgraph
			for(TargetNode nei: tmpgraphmaps.get(attackedtarget).getNeighbors())
			{
				if(tmpgraphmaps.keySet().contains(nei.getTargetid()))
				{
					System.out.println("targetid "+nei.getTargetid() + " , dist "+ tmpgraphmaps.get(attackedtarget).getDistance(nei));
					
					
					if(tmpgraphmaps.get(attackedtarget).getDistance(nei) <=5 && (nei.getTargetid() != 0))
					{
						attackcluster.add(nei.getTargetid());
					}
				}
				
			}
		}
		else
		{
			System.out.println("Attack target "+ attackedtarget +" is NOT in tmpraph ");
		}
		System.out.println();
		
		return attackcluster;
	
		
	}


	public static void dOWithAttackClusterTest(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int RADIUS, int slavelimit, int pathlimit) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int totalslaveiter = 0;
		long sumclustertime = 0;
		
		HashMap<Integer, Integer> clusterhistogram = new HashMap<Integer, Integer>();


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			/*if(iter==4)
			{
				System.out.println("xx");
			}
			*/
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			
			
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			
			double[] res = dOWithAttackCluster(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps, RADIUS, clusterhistogram, slavelimit, pathlimit);
			
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			totalslaveiter += res[6];
			sumclustertime += res[7];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//SecurityGameContraction.writeRes("DOClus", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);
		
		//writeClusterHist(clusterhistogram, ITER, nTargets);

		SecurityGameContraction.writeInFile("dOWithAttackCluster",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER, 
				totaltime/ITER, nTargets, totalslaveiter/ITER, sumclustertime/ITER, slavelimit, pathlimit);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	public static void dOWithAttackClusterTest2(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int RADIUS, int slavelimit, int pathlimit) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int totalslaveiter = 0;
		long sumclustertime = 0;
		
		HashMap<Integer, Integer> clusterhistogram = new HashMap<Integer, Integer>();


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			/*if(iter==4)
			{
				System.out.println("xx");
			}
			*/
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			
			
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			
			double[] res = dOWithAttackCluster2(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps, RADIUS, clusterhistogram, slavelimit, pathlimit);
			
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			totalslaveiter += res[6];
			sumclustertime += res[7];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//SecurityGameContraction.writeRes("DOClus", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);
		
		//writeClusterHist(clusterhistogram, ITER, nTargets);

		SecurityGameContraction.writeInFile("dOWithAttackCluster2",(int)sumfinaltargetsize/ITER, sumsol/ITER, 
				sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER, totaltime/ITER, nTargets, totalslaveiter/ITER, sumclustertime/ITER, slavelimit, pathlimit);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	public static void dOWithAttackClusterTest3(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int RADIUS, int slavelimit, int pathlimit) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int totalslaveiter = 0;
		long sumclustertime = 0;
		
		HashMap<Integer, Integer> clusterhistogram = new HashMap<Integer, Integer>();


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			/*if(iter==4)
			{
				System.out.println("xx");
			}
			*/
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			
			
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			
			double[] res = dOWithAttackCluster3(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps, RADIUS, clusterhistogram, slavelimit, pathlimit);
			
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			totalslaveiter += res[6];
			sumclustertime += res[7];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//SecurityGameContraction.writeRes("DOClus", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);
		
		//writeClusterHist(clusterhistogram, ITER, nTargets);

		SecurityGameContraction.writeInFile("dOWithAttackCluster3",(int)sumfinaltargetsize/ITER, sumsol/ITER,
				sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER, totaltime/ITER, nTargets, totalslaveiter/ITER, sumclustertime/ITER, slavelimit, pathlimit);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	public static void dOWithAttackClusterTest4(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int RADIUS, int slavelimit, int pathlimit) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int totalslaveiter = 0;
		long sumclustertime = 0;
		
		HashMap<Integer, Integer> clusterhistogram = new HashMap<Integer, Integer>();


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			/*if(iter==4)
			{
				System.out.println("xx");
			}
			*/
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			
			
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			
			double[] res = dOWithAttackCluster4(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps, RADIUS, clusterhistogram, slavelimit, pathlimit);
			
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			totalslaveiter += res[6];
			sumclustertime += res[7];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//SecurityGameContraction.writeRes("DOClus", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);
		
		//writeClusterHist(clusterhistogram, ITER, nTargets);

		SecurityGameContraction.writeInFile("dOWithAttackCluster4",(int)sumfinaltargetsize/ITER, sumsol/ITER,
				sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER, totaltime/ITER, nTargets, totalslaveiter/ITER, sumclustertime/ITER, slavelimit, pathlimit);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	
	
	public static void DOWithPACMANClusteringTest(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int RADIUS, int slavelimit, int pathlimit) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int totalslaveiter = 0;
		long sumclustertime = 0;
		
		HashMap<Integer, Integer> clusterhistogram = new HashMap<Integer, Integer>();


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			/*if(iter==4)
			{
				System.out.println("xx");
			}
			*/
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			
			
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			
			double[] res = DOWithPACMANClus(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps, RADIUS, clusterhistogram, slavelimit, pathlimit);
			
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			totalslaveiter += res[6];
			sumclustertime += res[7];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

		//	SecurityGameContraction.writeRes("DOWithPACMANClusteing", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);
		
		//writeClusterHist(clusterhistogram, ITER, nTargets);

		SecurityGameContraction.writeInFile("DOWithPACMANClusteing",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER,
				sumsolvtime/ITER, sumslavetime/ITER, totaltime/ITER, nTargets, totalslaveiter/ITER, sumclustertime/ITER, slavelimit, pathlimit);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	public static void DOWithSplitPACMANClusteringTest(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int RADIUS, int slavelimit, int pathlimit) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int totalslaveiter = 0;
		long sumclustertime = 0;
		
		HashMap<Integer, Integer> clusterhistogram = new HashMap<Integer, Integer>();


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			/*if(iter==4)
			{
				System.out.println("xx");
			}
			*/
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			
			
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = DOWithSplitPACMANClus(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps, RADIUS, clusterhistogram, slavelimit, pathlimit);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			totalslaveiter += res[6];
			sumclustertime += res[7];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//SecurityGameContraction.writeRes("DOWithSplitPACMANClusteing", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER, slavelimit, pathlimit);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);
		
		//writeClusterHist(clusterhistogram, ITER, nTargets);

		SecurityGameContraction.writeInFile("DOWithSplitPACMANClusteing",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER,totaltime/ITER,
				nTargets, totalslaveiter/ITER, sumclustertime/ITER, slavelimit, pathlimit);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	public static void writeClusterHist(HashMap<Integer,Integer> clusterhistogram, int ITER, int nTargets) 
	{


		try
		{
			
			File f = new File("result\\cluster-hist"+nTargets+".csv");
			 
			 if(f.exists())
			 {
				 f.delete();
				 f.createNewFile();
			 }
			
			
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("result\\cluster-hist"+nTargets+".csv"),true));
			
			for(Integer stid: clusterhistogram.keySet())
			{
				double size = clusterhistogram.get(stid)/ITER;
				pw.append(stid+","+size+"\n");
			}
			
			
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			//pw.append(expno+","+nTargets+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," +slavetime+","+ totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}





	}


	public static void naiveClusetringTest() {
		 
		
		int N = 8;
		int k = 2;
		
		int nTargets = N*N;
		
		int blocksize = k*k;
		
		
		int ncluster = nTargets/blocksize;
		
		
		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[ncluster];
		
		for(int i=0; i<ncluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}
		
		
		for(int c=0; c<ncluster; c++)
		{
			
			if(c==8)
			{
				System.out.println("x");
			}
			
			int starttarget = 0;
			
			if(c==0)
			{
				starttarget = 0;
			}
			else if(c%blocksize==0)
			{
				starttarget = c*blocksize;
			}
			else
			{
				// get the k-1th targetid for previous cluster c-1
				int id = clusters[c-1].get(k-1);
				starttarget = id+1;
			}
			
			
			int row = starttarget/N;
			int col = starttarget % N;
			
			
			for(int i=row; i<(row+k); i++)
			{
				for(int j=col; j<(col+k); j++)
				{
					int id = i*N + j;
					clusters[c].add(id);
				}
			}
			
			
			
		}
		
		printClusters(clusters);
		
		
		
	}


	public static String format(double avgsol) {
		// TODO Auto-generated method stub
		return new DecimalFormat("#00.000").format(avgsol);
	}

	
	
	

}



