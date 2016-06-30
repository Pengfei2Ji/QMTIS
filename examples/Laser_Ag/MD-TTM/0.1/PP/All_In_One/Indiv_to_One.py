#!/usr/bin/python

from itertools import izip

One_Output_File = open ('Gnuplot_Output.dat', 'w')
#Dens_Output_File = open ('Gnuplot_Dens_Output.dat', 'w')
#PressTens_Output_File = open ('Gnuplot_PressTens_x_Output.dat', 'w')
#Temp_Output_File = open ('../Gnuplot_Temp_Output.dat', 'w')


for j in range (0, 400):

#	IMD_Den_Indiv = open("../../Ag_Output.out.%d.dens" % j, "rU")
#	IMD_Str_Indiv = open("../../Ag_Output.out.%d.presstens" % j, "rU")

	with  open("../../Ag_Output.out.%d.dens" % j, "rU") as Dens_Indiv, open("../../Ag_Output.out.%d.presstens" % j, "rU") as PressTens_Indiv, open("./TTM_320_Round/Ag_Output.320.out.%d.ttm" % j, "rU") as Temp_Indiv:

		ii = 1
	
		for Line_Dens, Line_Presstens, Line_Temp in izip(Dens_Indiv, PressTens_Indiv, Temp_Indiv):
			if ii > 7:
				Col_Dens = Line_Dens.split()
				Col_PressTens = Line_Presstens.split()
				Col_Temp = Line_Temp.split()		

				Dens = float(Col_Dens[3])*1.661*107.8682/10.49
				PressTens = float(Col_PressTens[3])*160
				Te = float(Col_Temp[1])
				Tl = float(Col_Temp[2])
		
				if Dens < 0.1:
					PressTens = 0
					Te = 0
					Tl = 0 
						
#				Dens_Output_File.write('%d    %d    %13.8f\n' % ((j+1), (ii-7), Dens))
#				PressTens_Output_File.write('%d    %d    %13.8f\n' % ((j+1), (ii-7), PressTens))
				One_Output_File.write('%d    %d    %13.8f    %13.8f    %13.8f    %13.8f\n' % ((j+1), (ii-7), Dens, PressTens, Te, Tl))

			ii += 1
#		Dens_Output_File.write('\n')
#		PressTens_Output_File.write('\n')
		One_Output_File.write('\n')

	j += 1

#Dens_Output_File.close()
#PressTens_Output_File.close()

One_Output_File.close()
