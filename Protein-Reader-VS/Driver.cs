// ---------------------------------------------------------------------------
// File:		Reader.cs
// Project:		Analyzing SARS-COV2 spike proteins
// Author:		Jeffrey Richards, jmidrichards@gmail.com
// Creation:	05/19/2021
// ---------------------------------------------------------------------------
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using XPlot.Plotly;
using CsvHelper;
using System.Globalization;

namespace SpikeProteinProject
{
    class Driver
    {

        public static void Main(string[] args)
        {
			string[] proteins = Directory.GetFiles(@"..\..\..\iCompBio-Summer2021\Proteins-PDB", "*.*", SearchOption.AllDirectories);
				
			foreach (string proteinPath in proteins)
			{
				string proteinName = "6xra";
				Reader spike = new Reader($@"..\..\..\iCompBio-Summer2021\Proteins-PDB\{proteinName}.pdb");

				spike.ReadAtoms();
				var coordList = spike.coords.ToList<Coordinate>();
				coordList.RemoveAll(item => item == null);

				//Write coord arrays into a csv so that we can read in Python
				WriteCoords($@"..\..\..\iCompBio-Summer2021\Coordinates\{proteinName}.csv", coordList);

				Console.WriteLine(coordList.Count);
			}
			Console.ReadLine();
			Console.ReadLine();

            //string input = "";
            //while (input != "exit" && input != "show")
            //{
            //    Console.WriteLine("Enter 'show' to generate the plot for this protein, or 'exit' to cancel.");
            //    input = Console.ReadLine().ToLower();
            //}
            //if (input == "show")
            //    ShowScatter3d(spike);

        }

        public static void WriteCoords(string path, List<Coordinate> coords)
		{
            using (StreamWriter file = new StreamWriter(path))
            using (var csv = new CsvWriter(file, CultureInfo.InvariantCulture))
			{
                csv.WriteRecords(coords.ToList<Coordinate>());
			}
        }

		public static void ShowScatter3d(Reader reader)
		{
			ShowScatter3d(new Reader[] { reader });
		}

		public static void ShowScatter3d(Reader[] readers, string[] colors = null)
		{
			if (colors == null)
				colors = new string[] { "#1f77b4", "#9467bd", "#bcbd22", "#FD1E23", "#B46700", "#E7427C" };
			int maxNumPlots = 6;

			//Set up basic Scatter3d graph layout
			var plots = new Scatter3d[readers.Length];
			var labels = new string[] { "x", "y", "x" };
			var layout = new Layout.Layout();
			layout.title = "Proteins -";
			layout.autosize = false;
			layout.margin = new Margin();
			layout.margin.t = 60;

			for (int i = 0; i < readers.Length && i < maxNumPlots; i++)
			{
				Scatter3d plot = new Scatter3d();
				plot.x = TakeCol(readers[i].coords, "X");
				plot.y = TakeCol(readers[i].coords, "Y");
				plot.z = TakeCol(readers[i].coords, "Z");
				plot.mode = "lines";

				plot.line = new Line();
				plot.line.color = colors[i % colors.Length];
				plot.line.width = 4;

				layout.title += " " + readers[i].Name;
				plots[i] = plot;                //add current plot to the graph
			}

			//Create the graph using all input plots
			PlotlyChart chart = new PlotlyChart();
			chart.Plot(plots, layout, labels);
			chart.WithWidth(700);
			chart.WithHeight(500);
			chart.Show();

			//gives 275 KB html file containing 3d plot of the given proteins
		}

		public static float[] TakeCol(Coordinate[] coords, string col)
		{
			float[] retArray = new float[coords.Length];
			switch (col)
			{
				case "X":
					for (int i = 0; i < coords.Length; i++)
						retArray[i] = coords[i].X;
					break;
				case "Y":
					for (int i = 0; i < coords.Length; i++)
						retArray[i] = coords[i].Y;
					break;
				case "Z":
					for (int i = 0; i < coords.Length; i++)
						retArray[i] = coords[i].Z;
					break;
			}
			return retArray;
		}

	}
}
