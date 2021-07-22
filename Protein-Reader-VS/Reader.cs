// ---------------------------------------------------------------------------
// File:		Reader.cs
// Project:		Analyzing SARS-COV2 spike proteins
// Author:		Jeffrey Richards, jmidrichards@gmail.com
// Creation:	05/19/2021
// ---------------------------------------------------------------------------
using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpikeProteinProject
{
	class Reader
	{
		public string Name { get; private set; }
		private string[] Lines { get; set; }
		public Coordinate[] coords { get; set; }

		public Reader(string path)
		{
			Lines = File.ReadAllLines(path);
			Console.WriteLine(Lines[0]);
			Name = path.Substring(path.Length - 8, 4).ToUpper();

			coords = new Coordinate[Lines.Length / 5];
		}

		public void ReadAtoms(string atom = "A")
		{
			atom = atom.ToUpper();

			for (int i = 0, j = 0; i < Lines.Length; i++)   //iterate through every line of the file
			{
				string[] line = Lines[i].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

				if (line[0] == "ATOM")
				{
					if (line[4].ToCharArray().Length > 1)   //if col 4/5 are 'A1000' instead of 'A 999'
					{
						if (line[4].ToCharArray()[0].ToString() == atom && line[2] == "CA")
						{
							coords[j] = new Coordinate(Single.Parse(line[5]), Single.Parse(line[6]), Single.Parse(line[7]));
							j++;                            //How many coordinates we've put into our arrays so far
						}
					}
					else if (line[4] == atom && line[2] == "CA")
					{
						coords[j] = new Coordinate(Single.Parse(line[6]), Single.Parse(line[7]), Single.Parse(line[8]));
						j++;
					}

				}//end if ATOM

			}//end for every line

		}

		public void ReadAtoms(int atom)
		{

		}

		public bool WriteCoords(string path)
		{
			return false;
		}

	}
}
