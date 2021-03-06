% INPUT:
% -filename is the file that contains the instances generated by netgen
% -mode is the execution mode:
%                             -'precond' compares the peconditioned system and the original one
%                             -'minres' compares CustomGMRES and the MATLAB function minres
%                             -'residual' plots the relative residual wrt the number of iterations
%                             -'rate' plots the rate of convergence
% -generate is a boolean variable which indicates if the function has to generate or to load the vector D                            
% -distribution is the probability distribution of data in the vector D choosen between 7 fixed distributions

function main(filename, mode, generate, distribution)

    if strcmp(mode, 'precond')
        precond = true;
    else
        precond = false;
    end
    
    prompt = "Insert the number of nodes of the graph: ";
    nodi = input(prompt);
    
    flag_nop = true; %flag to stop the algorithm (no precond.) when it reaches the delta
    flag_p = false; %%flag to stop the algorithm (preconditioned) when it reaches the delta
    a = 0; %perturbation S+a*I

    if ~isnumeric(nodi)
        ME = MException('%s is not a number',nodi);
        throw(ME)
    end

    [E, b, c] = netgenreader(filename, nodi); %read the data
    e = size(E,2); % #edges
    
    b=-b;
    b_tilde = [b; c]; % given b and c we build b tilde
    
    b_norm = norm(b_tilde);
    n = nodi + e;
    prompt = "Insert delta value. Insert 0 for the default value: ";
    in = input(prompt);

    if in == 0
        delta = 10^(-10);
    else 
        if isnumeric(in) && in > 0
            delta = in;
        else
            ME = MException('%s is not a number or is negative',in);
            throw(ME)
        end
    end

    if generate %chech if we have to generate D or to load it
        rng(46273); %seed
        switch distribution
        case 'Gamma'
            D = random('Gamma', 5,1, e, 1);
            save('D_gamma', 'D');
        case 'Beta75'
            D = random('Beta', 0.75, 0.75, e, 1);
            save('D_beta75', 'D');
        case 'Beta44'
            D = random('Beta', 4, 4, e, 1);
            save('D_beta44', 'D');
        case 'Chi'
            D = random('chi2', 4, e, 1);
            save('D_chi2', 'D');
        case 'MixBin'
             D=random('Binomial', 20, 0.7, e,1)+9;
             perm = randperm(e, round(e*0.7));
             D(perm) = random('Binomial', 20, 0.7, round(e*0.7),1);
             save('D_mixbin', 'D');
        case 'Uniform'
            D = rand(e,1);
            save('D_uniform', 'D');
        case 'Ill'
             D=rand(e,1);
             perm = randperm(e, round(e*0.4));
             D(perm)= D(perm)*10^-5;
             save('D_ill', 'D');
        end
    else
        load(distribution)
    end
    

    %INITIALIZING THE VARIABLES FOR EACH MODE
    switch mode 
        case 'minres'
            % Building A in sparse format
            A = sparse([1:e], [1:e], D, nodi+e, nodi+e);
            A(e+1:end,1:e) = E;
            A(1:e, e+1:end) = E';
            jump = round(sqrt(n)/5); %fix the iterations jump
        case 'precond'
            jump = round(sqrt(n)/5);
        otherwise
            jump = 1;
    end
    
    num_iterations = [50:jump:round(8*sqrt(n))]; %fix the iterations we want to do

    %initializing variables
    x_custom = zeros(length(b_tilde),length(num_iterations)); % solution 
    % of customized GMRES with different number of iterations 
    t_nop = zeros(length(num_iterations),1); % time spent by customized 
    % GMRES with different number of iterations
    residual_nop = zeros(length(num_iterations), 1);
    
    if strcmp(mode, 'rate')
        rate = zeros(length(num_iterations)-1, 1);
    end
    
    q_nop = zeros(length(num_iterations)-1, 1);
    
    if precond
        flag_p = true;
        prompt="Insert the value of perturbation of S: a = ";
        a = input(prompt);
        residual_p = zeros(length(num_iterations),1);
        residual_px = zeros(length(num_iterations),1);

        t_p2 = zeros(length(num_iterations),1);
        t_px = zeros(length(num_iterations),1);
        q_p = zeros(length(num_iterations),1);
    end
    
    %START THE ACTUAL EXECUTION
    for i = 1:length(num_iterations)
        fprintf('Iteration number %d (m=%d)\n', i, num_iterations(i));
        
         [x_custom(:, i), q_nop(i),q_p(i), residual_nop(i),residual_p(i), t_nop(i), t_p2(i), t_px(i)] = find_sol(D, E, b, c, num_iterations(i), flag_p, flag_nop, a); 
           
         if flag_p %computing the Preconditioned Original residual and time
            tic;
            Ax = [D.*x_custom(1:length(D),i) + E'*x_custom(length(D)+1:end,i); E*x_custom(1:length(D),i)];
        
            residual_px(i) = norm(Ax-b_tilde)/b_norm; % compute the residual
            
            t_temp = toc;
            t_px(i) = t_px(i) + t_temp;
         end

        if strcmp(mode, 'rate') && i > 1
            rate(i-1) = residual_nop(i)/residual_nop(i-1);
        end
 
        if i > 1 %check the stopping criterion
            del_nop = abs(residual_nop(i)-residual_nop(i-1));
            if del_nop < delta && flag_nop
                disp("CustomGMRES stopped since the delta was close to 0 (not preconditioned)");
                flag_nop = false;
                position_nop = i;
            end
            
            if strcmp(mode, 'precond')
                del_p = abs(residual_p(i)-residual_p(i-1));
                if del_p < delta && flag_p
                    disp("CustomGMRES stopped since the delta was close to 0 (preconditioned)");
                    flag_p = false;
                    position_p = i;
                end
            end

        end  

        %if both precond. and not precond. stopped (or not precond stopped
        %and mode is not 'precond') we stop
        if ~flag_nop && ~flag_p  
            break;
        end
    end 
    
    if flag_nop
        position_nop = i;
    end
   
    switch mode %plots
        case 'precond'
            if flag_p
                position_p = i;
            end
            scatter(q_nop(1:position_nop),t_nop(1:position_nop),80,"red", "filled");
            hold on;
            scatter(q_p(1:position_p),t_p2(1:position_p), 80,"magenta","filled");

            h(1) = scatter(residual_nop(1:position_nop), t_nop(1:position_nop), 35, 'blu','square', 'filled'); % no precond
            h(2) = scatter(residual_p(1:position_p), t_p2(1:position_p), 35, 'green', 'filled'); % precond transformed
            h(3) = scatter(residual_px(1:position_p),t_px(1:position_p), 45, [0.4940 0.1840 0.5560],'filled'); % precond original
            legend(h, ["Non Preconditioned","Preconditioned Transformed","Preconditioned Original"]);
            legend('Location', 'northwest');

            xlabel('Relative Residual');
            ylabel('Time (s)');

            set(gca, 'XScale', 'log');
        case 'rate'
            ratefunc(num_iterations(1:i-1),rate(1:i-1))
        case 'residual'
            residualfunc(num_iterations(1:i), residual_nop(1:i))
        case 'minres'
            c = 1:i;
            h(1)=scatter(q_nop(1:i),t_nop(1:i),80,"red", "filled");
            hold on;
            h(2)=scatter(residual_nop(1:i), t_nop(1:i), 25, c, 'filled');
            colormap(winter);
            r_minres = zeros(i,1);
            t_minres = zeros(i,1);
            
            for j = 1:i
                tic;
                [~,~,r_minres(j)] = minres(A, b_tilde, 1e-15, num_iterations(j));
                t_minres(j) = toc;
            end
            
            h(3)=scatter(r_minres, t_minres, 60, [0.3 0.2 0.6], 'square', 'filled');

            legend(h, ["Lower Bound","CustomGMRES","MINRES"]);
            set(gca, 'XScale', 'log');       
    end
end
