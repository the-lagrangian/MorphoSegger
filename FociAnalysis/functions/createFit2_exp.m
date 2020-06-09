function [fitresult, g_fit, t_t] = createFit2_exp(x, a,p_fit,exp_cut)
%CREATEFIT1(X,A)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: a
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 30-Mar-2020 19:27:23


%% Fit: 'untitled fit 1'.

%x=tiempo;
%a=l_cell;

[~, loc_p]=findpeaks(a);

        if isempty(loc_p)==1

            xx=x;

            aa=a;

        else
            
            xt1=x(1:loc_p(1));
            aa_t1=a(1:loc_p(1));
            
            xt2=x(loc_p(1)+1:end);
            aa_t2=a(loc_p(1)+1:end);
            
            
                if length(xt1)>length(xt2)
                
                    xx=xt1;
                    aa=aa_t1;
                
                else
                    xx=xt2;
                    aa=aa_t2;
            
                end    
      
                
        end

loc_a=find(xx>exp_cut); % Takes only 65 pixels onwards to fit the exponential

xx=xx(loc_a);
aa=aa(loc_a);      
        
t_t=0;

if length(xx)<p_fit
    
    t_t=1;
end

[xData, yData] = prepareCurveData( xx, aa );



% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Display = 'Off';
%opts.StartPoint = [34.182834814215 0.00653034302260214];

% Fit model to data.
[fitresult, g_fit] = fit( xData, yData, ft, opts );

%[xData, yData] = prepareCurveData( x, a );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'a vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'x', 'Interpreter', 'none' );
% ylabel( 'a', 'Interpreter', 'none' );
% grid on

 if g_fit.rsquare<0.95
    
    t_t=1;
 end
 

