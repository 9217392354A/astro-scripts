#program to compare values 

import numpy as np



rho_g_s_coma = [0.0058083497, 0.0055145699]
rho_d_s_coma = [9.90384e-05, 7.04471e-05]

rho_g_s_filament = [0.040263, 0.037849]
rho_d_s_filament = [0.000596, 0.000439]

rho_g_s_virgo = [0.040263, 0.037849]
rho_d_s_virgo = [0.000596, 0.000439]

rho_g_s_field = [0.3425459012, 0.0936941533]
rho_d_s_field = [0.0010307244, 0.00334188]


names = ['coma', 'filament', 'virgo', 'field']
gas_to_stars  = [rho_g_s_coma, rho_g_s_filament, rho_g_s_virgo, rho_g_s_field]
dust_to_stars = [rho_d_s_coma, rho_d_s_filament, rho_d_s_virgo, rho_d_s_field]


if True:
	l = []
	print names
	#gas grid
	for j in range(4):
		y_name = names[j]
		y = gas_to_stars[j][0]
		dy = gas_to_stars[j][1]

		line = []

		for i in range(4):
			x_name = names[i]
			x = gas_to_stars[i][0]
			dx = gas_to_stars[i][1]

			if x_name == y_name:
				diff = '-'

			else:
				diff = np.round(abs(x-y)/np.sqrt(dx**2 + dy**2), decimals=0)

			line.append(diff) #y_name, 'vs', x_name, diff


		print y_name, line


if True:
	l = []
	print names
	#gas grid
	for j in range(4):
		y_name = names[j]
		y = dust_to_stars[j][0]
		dy = dust_to_stars[j][1]

		line = []

		for i in range(4):
			x_name = names[i]
			x = dust_to_stars[i][0]
			dx = dust_to_stars[i][1]

			if x_name == y_name:
				diff = '-'

			else:
				diff = np.round(abs(x-y)/np.sqrt(dx**2 + dy**2), decimals=0)

			line.append(diff) #y_name, 'vs', x_name, diff


		print y_name, line