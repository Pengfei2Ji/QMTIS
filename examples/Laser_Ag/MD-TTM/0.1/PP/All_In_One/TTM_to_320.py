#!/bin/python

for j in range (0, 400):

	IMD_TTM_SP_Au = open('../../Ag_Output.out.%d.ttm' % j, 'rU')
	TTM_to_320_Output_File = open ('./TTM_320_Round/Ag_Output.320.out.%d.ttm' %j, 'w')

	ii = 1
	
	Line_Num_Now = 0.0
	Line_Num_Pre = 0.0

	Te_Pre = 0.0
	Tl_Pre = 0.0

	Te_Add = 0.0
	Tl_Add = 0.0

	for line in IMD_TTM_SP_Au:
		if ii == 1:
			TTM_to_320_Output_File.write('# Line_001 \n')
			TTM_to_320_Output_File.write('# Line_002 \n')
			TTM_to_320_Output_File.write('# Line_003 \n')
			TTM_to_320_Output_File.write('# Line_004 \n')
			TTM_to_320_Output_File.write('# Line_005 \n')
			TTM_to_320_Output_File.write('# Line_006 \n')
			TTM_to_320_Output_File.write('# Line_007 \n')

		else:
			Col_TTM_SP_i = line.split()
			Te = float(Col_TTM_SP_i[4])*11605
			Tl = float(Col_TTM_SP_i[5])*11605
			
			Line_Num_Now = (round( (ii-1)/288.0*320 ) - 1)

			if ( (Line_Num_Now - Line_Num_Pre) > 1):
				Line_Num_Add = Line_Num_Pre +1
				Te_Add = (Te_Pre + Te)/2.0
				Tl_Add = (Tl_Pre + Tl)/2.0
				TTM_to_320_Output_File.write('%13.8f    %13.8f    %13.8f\n' % (Line_Num_Add, Te_Add, Tl_Add))

			TTM_to_320_Output_File.write('%13.8f    %13.8f    %13.8f\n' % (Line_Num_Now, Te, Tl))

			Te_Pre = Te
			Tl_Pre = Tl
			Line_Num_Pre = Line_Num_Now

#			if Tl == 0:
#				Tl = Te

#			Te_Output_File.write('%d    %d    %s\n' % ((j+1), (ii-1), Te))
#			Tl_Output_File.write('%d    %d    %s\n' % ((j+1), (ii-1), Tl))
#			Tl_Melt_Output_File.write('%d    %d    %s\n' % ((j+1), (ii-1), 1234.93))
#			Tl_Output_File.write('%d    %d    %s\n' % ((j+1), (ii-1), 2435))


		ii += 1
#	Te_Output_File.write('\n')
#	Tl_Output_File.write('\n')

	TTM_to_320_Output_File.close()

	j += 1

#Te_Output_File.close()
#l_Output_File.close()
#Tl_Melt_Output_File.close()
#Tl_Evap_Output_File.close()
