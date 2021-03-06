package groupingtargets;

import java.awt.Color;
import java.awt.GridLayout;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JButton;
import javax.swing.JFrame;

import cs.Interval.Abstraction.SecurityGameAbstraction;
import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;

class ButtonGrid {





	JFrame frame=new JFrame(); //creates frame
	JButton[][] grid; //names the grid of buttons
	
	JFrame clusterframe=new JFrame(); //creates frame
	JButton[][] clustergrid; //names the grid of buttons
	
	
	HashMap<Integer, SuperTarget> sts;
	HashMap<Integer, TargetNode> targetmaps;
	ArrayList<Color> colors = new ArrayList<Color>();

	public ButtonGrid( HashMap<Integer, TargetNode> targetmaps, HashMap<Integer, SuperTarget> sts, String title)
	{ //constructor


		this.sts = sts;
		this.targetmaps = targetmaps;
		
		this.clusterframe.setTitle(title);
		
		
		
		
	}
	
	
	public void drawPayoffGrid(int nrow, int ncol)
	{
		frame.setLayout(new GridLayout(nrow,ncol)); //set layout
		grid=new JButton[nrow][ncol]; //allocate the size of grid
		

		//int c = 255;
		
		int targetid = 0;

		for(int y=0; y<ncol; y++)
		{
			for(int x=0; x<nrow; x++)
			{
				int au = (int)targetmaps.get(targetid).attackerreward;
				grid[x][y]=new JButton(String.valueOf(targetid)+"("+au+")"); //creates new button 
				grid[x][y].setOpaque(true);
				grid[x][y].setBorderPainted(false);


				
				
				au *= 10;
				
				

				grid[x][y].setBackground(new java.awt.Color(155+au,155-au,155-au));
				frame.add(grid[x][y]); //adds button to grid
				
				targetid++;
				
			}
		}
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack(); //sets appropriate size for frame
		frame.setVisible(true); //makes frame visible
	}
	
	public void drawCluster(int nrow, int ncol)
	{
		
		
		HashMap<Integer, java.awt.Color> clustercolor = new HashMap<Integer, java.awt.Color>();
		
		clusterframe.setLayout(new GridLayout(nrow,ncol)); //set layout
		clustergrid=new JButton[nrow][ncol]; //allocate the size of grid

		//int c = 255;
		
		int targetid = 0;

		for(int y=0; y<ncol; y++)
		{
			for(int x=0; x<nrow; x++)
			{
				
				
				int clusid = getClusId(targetid);
				

				//if(clusid!=-1)
				//{
					int au = (int)targetmaps.get(targetid).attackerreward;

					
					clustergrid[x][y]=new JButton(String.valueOf(targetid+"("+clusid+")")); //creates new button 
					clustergrid[x][y].setOpaque(true);
					clustergrid[x][y].setBorderPainted(false);
					//clustergrid[x][y].
					


					//int au = (int)targetmaps.get(targetid).attackerreward;

					au *= 10;

					java.awt.Color clusc = Color.WHITE;


					// if cluster has more than one node than assign color
					if(clusid>=0 && sts.get(clusid).nodes.size()>1)
					{
						//check if it already has a color assignment
						if(clustercolor.keySet().contains(clusid))
						{
							clusc = clustercolor.get(clusid);
						}
						else
						{
							// get a random color 
							int r = SecurityGameContraction.randInt(0,255);
							int g = SecurityGameContraction.randInt(0, 255);
							int b = SecurityGameContraction.randInt(0, 255);

							Color tmp = new java.awt.Color(r,g,b);
							clustercolor.put(clusid, tmp);
							clusc = tmp;
						}
						//clustergrid[x][y].setBackground(clusc);
						
						
						clustergrid[x][y].setForeground(clusc);
						clustergrid[x][y].setBackground(new java.awt.Color(155+au,155-au,155-au));
					}
					else if(clusid>=0 && (sts.get(clusid).nodes.size()==1))
					{
						clustergrid[x][y].setBackground(new java.awt.Color(155+au,155-au,155-au));
					}
					else
					{
						clustergrid[x][y].setBackground(Color.BLACK);
						clustergrid[x][y].setForeground(Color.WHITE);
					}

					// if not then assign white

					//au*=10;
					clustergrid[x][y].setOpaque(true);
					clusterframe.add(clustergrid[x][y]); //adds button to grid
					clusterframe.setOpacity(1);
				//}
				
				
				targetid++;
				
			}
		}
		clusterframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		clusterframe.pack(); //sets appropriate size for frame
		clusterframe.setVisible(true); //makes frame visible
	}
	/* public static void main(String[] args) {
           new ButtonGrid(10,10);//makes new ButtonGrid with 2 parameters
   }
	 */


	private int getClusId(int targetid) {
		
		
		for(SuperTarget st: sts.values())
		{
			for(TargetNode t: st.nodes.values())
			{
				if(t.getTargetid() == targetid)
					return st.stid;
			}
		}
		
		
		return -1;
	}

}
