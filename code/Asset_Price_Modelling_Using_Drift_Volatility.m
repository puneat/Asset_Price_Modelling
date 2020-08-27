clc;clear all;close all
%%
%{
 Importing data in the form of .CSV and storing as table in MATLAB
%}
T_msft = readtable('C:\Users\punit\Desktop\OMEGA\MSFT_daily.csv','PreserveVariableNames',1);
T_aapl = readtable('C:\Users\punit\Desktop\OMEGA\SBICARD_daily.csv','PreserveVariableNames',1);
test = readtable('C:\Users\punit\Desktop\OMEGA\MSFT_test.csv',1);
%%
%{
Using a predefined function all_to_close (available at the bottom of the script)
to extract the closing price of all the days from the loaded data file
%}
msft_close=all_to_close(T_msft,4);
aapl_close=all_to_close(T_aapl,4);
test_close=all_to_close(test,4);
%%
%{
 Calculating and storing the parameters of the daily closing price, namely,
1. Mean
2. Standard Devaiation
3. Variance
4. Drift

Further, the daily returns are also stored
%}
[msft_daily,msft_params]=calc_parameters(msft_close);
[aapl_daily,aapl_params]=calc_parameters(aapl_close);
%%
%{
 We define the number of simulations to run and number of days to predict
 in the future. Follwing which, we run the monte-carlo simulations
%}
days=39;
number_of_sims=1000;
msft_next_day_price=mcs_next_price(msft_close,msft_params,days,number_of_sims);
aapl_next_day_price=mcs_next_price(aapl_close,aapl_params,days,number_of_sims);
%%
%{
 Plotting the result of simulations 
%}
figure('Color','white','Name','Monte-Carlo Simulation','NumberTitle','off');
subplot(2,1,1)
plot(msft_next_day_price)
grid on
title('MSFT Monte Carlo Simulation')
xlabel('Days')
ylabel('Price')
% If we want to test against actual data
%hold on
%plot(test_close, 'LineWidth',3,'Color','black')
subplot(2,1,2)
plot(aapl_next_day_price)
grid on
title('SBI-CARD Monte Carlo Simulation')
xlabel('Days')
ylabel('Price')

%%
%{
 Plotting the histogram of the simulations and thus confirming if the
 resultant is normally distributed and finding the range of the most likely
 price estimates.
%}
figure('Color','white','Name','Frequency of Results','NumberTitle','off');
subplot(2,1,1)
hist(msft_next_day_price)
grid on
title('MSFT')
xlabel('Price')
ylabel('Frequency')
subplot(2,1,2)
hist(aapl_next_day_price)
grid on
title('SBI-CARD')
xlabel('Price')
ylabel('Frequency')
%%
%{
all_to_close(T,close_index)

DEF: Function to extract the closing price column from a table.

INPUT:
    T = Table of size (h,w) from which closing price is to be extracted.
    close_index = The column position of the closing price in the table.
 
OUTPUT:
    close_array = Array of size (h,1) containing the closing price vector.

%}
function [close_array] = all_to_close(T,close_index)
    T=T(:,2:6);
    T=table2array(T);
    close_array=T(:,close_index);
end
%%
%{
calc_parameters(close)

DEF: Function to calculate and return the various statistical parameters
of the given closing price vector.

INPUT:
    close = closing price vector of size (n,1).
 
OUTPUT:
    return_array = Array of size (4,1) containing the following statistical
                  parameters,
                  1. Mean of the vector
                  2. Standard Deviation of the vector
                  3. Variance of the vector
                  4. Drift of the vector
    daily_return = Array of size (n,1) containing the daily returns with
                   the n'th index as base.
%}

function [daily_return,return_array] = calc_parameters(close)
% Calculate the size of array (n,m)
    size_arr=size(close)
% Create an array for daily returns of size (n,1) and fill it with zeros.
    daily_return=zeros(size_arr(1),1);
% Create an array for statistical parameters of size (4,1) and fill it with zeros.
    return_array=zeros(4,1);
%{
    Calculate the daily returns for each day using the following formula,
    return_today= ln(price_today/price_yesterday)
%}
    for i=1:size_arr(1)-1
        daily_return(i)=log(close(i)/close(i+1));
    end
 % Variance
    var_return=var(daily_return);
    return_array(3)=var_return
 % Mean
    mean_return=mean(daily_return);
    return_array(1)=mean_return
 %{
    Calculate the overall drift using the following formula,
    drift= mean - variance/2
%}   
    drift=mean_return-(var_return/2);
    return_array(4)=drift
% Standard Deviation    
    std_return=std(daily_return);
    return_array(2)=std_return

end

%%
%{
mcs_next_price(close_price,params,days,number_of_sims)

DEF: Function to calculate and return the next day price of an equity
starting from a base point. Uses Monte-Carlo simulations to run the model n
number of times to cover all possible scenarios and incorporate randomness
of the given closing price vector.

INPUT:
    close_price = closing price vector of size (n,1).
    params = Array of size (4,1) containing the following statistical
                  parameters,
                  1. Mean of the vector,
                  2. Standard Deviation of the vector,
                  3. Variance of the vector,
                  4. Drift of the vector.

    days = Number of days in future to predict the price for

    number_of_sims = Number of times to run the simulation for
 
OUTPUT:
    next_day_price = Array of size (days,1) containing the future prices                 
%}
function next_day_price= mcs_next_price(close_price,params,days,number_of_sims)
% Create an array for next day prices of size (days,1) filled with zeros.
    next_day_price=zeros(days,number_of_sims);
% The base price for future predictions is chosen as the last price in the
% closing price vector.
    next_day_price(1,:)=close_price(1,:);
%{
The simulation is run for number_of_sims times and for the specified number
of days. Volatility of the asset is modelled by generating a random value
at each iteration and transforming the value through a normal cumalative 
probability curve. The next day price, is hence calculated as,
    price_next_day = prev_day_price * exp(drift+volatility)
%}
    for j=1:number_of_sims
        for i=2:days
            random_value=params(2) * norminv(rand());
            next_day_price(i,j)=next_day_price(i-1,j)*exp(params(4)+random_value);
        end
    end
end
%%
    
 
