%{ Usage: Mixture_DispStatus
  This will print a message to the command window with information about
  the status of the current Multivariate Mixture Modeler run. It replaces
  the nice status panel that was previously displayed in the GUI. The data,
  stored in the MMM object, that is printed includes:
  lblCurrRun = Current Replication # of Total Replications
  lblCurrIC = Either the name of the initialization function currently
     running, or the name of the information criterion currently running
  lblCurrk= Current k of Maximum K to fit
  lblTimRep = Processing time of the last completed replication
  lblTimK = Processing time of the last completed k run
  lblTimInit = Processing time of the last completed execution of the 
     initialization function
  lblTimIC = Processing time of the last execution of the optimization and
     information criterion function
  lblTimTot = Total elapsed processing time
  all processing times are in minutes

  Copyright (C) 2017 Prof. Hamparsum Bozdogan & J. Andrew Howe
%}

  
disp('*--- Status Update ---')
disp(['*    Current Replication: ',MMM.lblCurrRun])
disp(['*    Current k: ',MMM.lblCurrk])
disp(['*    Current IC: ',MMM.lblCurrIC])
disp('*Time of Last (m):')
disp(['*    ',MMM.lblTimRep])
disp(['*    ',MMM.lblTimK])
disp(['*    ',MMM.lblTimInit])
disp(['*    ',MMM.lblTimIC])
disp(['*Elapsed Time: ',MMM.lblTimTot])
disp('*---------------------')
pause(1)

%{
Copyright (C) 2017 Prof. Hamparsum Bozdogan & J. Andrew Howe

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}