#house buying

import numpy as np
import matplotlib.pyplot as plt


balance = 10**5
house_cost = 10**5
rent = 600.0*12.0

number_houses = np.floor(balance/ house_cost)


years = np.arange(2014,2014+40,1)

balances  = []
houses = []
rents = []
working_balance = []

net_worth = []

for year in years:
	working = balance
	if balance > house_cost: #buying possible
		#how many houses to buy
		n = np.int(np.floor(balance/house_cost))

		#add houses to current stock
		number_houses += n

		#remove money from balance for houses
		balance = balance - house_cost*n
	
	#get rent payed
	balance += rent*number_houses




	

	#save all
	balances.append(balance)
	working_balance.append(working)
	houses.append(number_houses)
	rents.append(rent*number_houses)
	net_worth.append(number_houses*house_cost+balance)
fig = plt.figure(figsize = (12.,12.),facecolor='w',edgecolor='w')
sub1 = fig.add_subplot(2,2,1)
sub2 = fig.add_subplot(2,2,2)
sub3 = fig.add_subplot(2,2,3)
sub4 = fig.add_subplot(2,2,4)

sub1.plot(years, np.array(balances, dtype = np.float)/10**6.)
sub1.plot(years, np.array(working_balance, dtype = np.float)/10**6.)
sub1.set_xlabel('Years')
sub1.set_ylabel('Balance (GBP millions)')

sub2.plot(years,houses)
sub2.set_xlabel('Years')
sub2.set_ylabel('Number of houses')

sub3.plot(years, np.array(rents, dtype = np.float)/10**6.)

sub3.set_xlabel('Years')
sub3.set_ylabel('Rental income per year (GBP millions)')

sub4.plot(years, np.array(net_worth, dtype = np.float)/10**6.)
sub4.set_xlabel('Years')
sub4.set_ylabel('Net worth (GBP millions)')

plt.show()

