function simlgca

    nodes_x=256;
    nodes_y=128;
    tsteps=400;
    
    nodes = zeros(nodes_x,nodes_y,6);
    
    v1 = [1; 0];
    v2 = [cos(pi/3); sin(pi/3)];
    v3 = [cos(2*pi/3); sin(2*pi/3)];
    v4 = [-1; 0];
    v5 = [cos(4*pi/3); sin(4*pi/3)];
    v6 = [cos(5*pi/3); sin(5*pi/3)];
    
    obstacle=zeros(nodes_x,nodes_y);
    
    for j=10:100
        obstacle(40,j)=1;
    end
    
    for i=1:nodes_x
        for j=(2:nodes_y - 1)
            if obstacle(i,j) ~=1
                current_cell=nodes(i,j,:);
                
                current_cell(1)=1;
                
                nodes(i,j,:)=current_cell;
            end
        end
    end    
                
    for t=1:tsteps
        for i=1:nodes_x
            for j=2:(nodes_y-1)
                if obstacle(i,j) ~=1
                    cell=nodes(i,j,:)
                    
                    num_particles=sum(cell);
                    
                    if (num_particles~=2) && num_particles~=3)
                        nodes(i,j,:)=cell
                        
                    elseif num_particles==3
                        if (cell(1)== cell(3))&&(cell(3)==cell(5))
                            nodes(i,j,:) =~ cell
                            
                        else
                            nodes(i,j,:) = cell
                            
                        end
                    else
                        p1 = find(cell, 1);
                        
                        if(p1>3)||(cell(p1 + 3) ~= 1)
                            nodes(i, j, :) = cell;
                            
                        else
                            random=rand
                            
                            if random < 0.5    
                                n_cell(1) = cell(6);
                                n_cell(2) = cell(1);
                                n_cell(3) = cell(2);
                                n_cell(4) = cell(3);
                                n_cell(5) = cell(4);
                                n_cell(6) = cell(5);

                            else             
                                n_cell(1) = cell(2);
                                n_cell(2) = cell(3);
                                n_cell(3) = cell(4);
                                n_cell(4) = cell(5);
                                n_cell(5) = cell(6);
                                n_cell(6) = cell(1);
                            end
                            
                            nodes(i, j, :) = n_cell;
                        end
                    end
                end
            end
        end
        
        for i = 1:nodes_x
            nodes(i, 1, :) = [nodes(i, 1, 4) nodes(i, 1, 5) nodes(i, 1, 6) nodes(i, 1, 1) nodes(i, 1, 2) nodes(i, 1, 3)];
            nodes(i, nodes_y, :) = [nodes(i, nodes_y, 4) nodes(i, nodes_y, 5) nodes(i, nodes_y, 6) nodes(i, nodes_y, 1) nodes(i, nodes_y, 2) nodes(i, nodes_y, 3)];
        end
        
        for i = 1:nodes_x
            for j = 1:nodes_y
                if (obstacle(i, j) == 1)
                    nodes(i, j, :) = [nodes(i, j, 4) nodes(i, j, 5) nodes(i, j, 6) nodes(i, j, 1) nodes(i, j, 2) nodes(i, j, 3)];
                end
            end
        end
        
        n_nodes = zeros(nodes_x, nodes_y, 6);
        
        for i = 1:nodes_x
            for j = 1:nodes_y
                
                cell = nodes(i, j, :);

              
                neighbor_x = 0;
                neighbor_y = 0;

                % Propagation in the 1-direction.
                neighbor_y = j;

                if i == nodes_x
                    neighbor_x = 1; 
                else
                    neighbor_x = i + 1;
                end

                n_cell = n_nodes(neighbor_x, neighbor_y, :);
                n_cell(1) = cell(1);
                n_nodes(neighbor_x, neighbor_y, :) = n_cell;

                % Propagation in the 2-direction.
                if j ~= nodes_y
                    neighbor_y = j + 1;

                    if mod(j, 2) == 0
                        if i == numnodes_x
                            neighbor_x = 1;
                        else
                            neighbor_x = i + 1;
                        end
                    else
                        neighbor_x = i;
                    end

                    n_cell = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell(2) = cell(2);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell;
                end

                % Propagation in the 3-direction.
                if j ~= nodes_y
                    neighbor_y = j + 1;

                    if mod(j, 2) == 1
                        if i == 1
                            neighbor_x = numnodes_x;
                        else
                            neighbor_x = i - 1;
                        end
                    else
                        neighbor_x = i;
                    end

                    n_cell = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell(3) = cell(3);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell;
                end

                % Propagation in the 4-direction.
                neighbor_y = j;

                if i == 1
                    neighbor_x = numnodes_x;
                else
                    neighbor_x = i - 1;
                end

                n_cell = n_nodes(neighbor_x, neighbor_y, :);
                n_cell(4) = cell(4);
                n_nodes(neighbor_x, neighbor_y, :) = n_cell;

                % Propagation in the 5-direction.
                if (j ~= 1)
                    neighbor_y = j - 1;

                    if (mod(j, 2) == 1)
                        if (i == 1)
                            neighbor_x = nodes_x;
                        else
                            neighbor_x = i - 1;
                        end
                    else
                        neighbor_x = i;
                    end

                    n_cell = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell(5) = cell(5);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell;
                end

                % Propagation in the 6-direction.
                if (j ~= 1)
                    neighbor_y = j - 1;

                    if (mod(j, 2) == 0)
                        if (i == nodes_x)
                            neighbor_x = 1;
                        else
                            neighbor_x = i + 1;
                        end
                    else
                        neighbor_x = i;
                    end

                    n_cell = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell(6) = cell(6);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell;
                end
            end
        end
        
        nodes = n_nodes;
        
    end
    
    
    grain_size = 8;
    grain_x = nodes_x / grain_size;
    grain_y = nodes_y / grain_size;
    
    av_vel_x_coords = zeros(1, grain_x * grain_y);
    av_vel_y_coords = zeros(1, grain_x * grain_y);
    av_vel_x_comps = zeros(1, grain_x * grain_y);
    av_vel_y_comps = zeros(1, grain_x * grain_y);
    
    currval = 1;
    for (i = 1:1:grain_x)
            % Calculate the lower and upper x-boundaries.
            x_bd_l = (i - 1)*grain_size + 1;
            x_bd_u = i*grain_size;
        for (j = 1:1:grain_y)
            % Calculate the lower and upper y-boundaries.
            y_bd_l = (j - 1)*grain_size + 1;
            y_bd_u = j*grain_size;      
            
            
            np = zeros(1, 6);
            np(1) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 1)));
            np(2) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 2)));
            np(3) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 3)));
            np(4) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 4)));
            np(5) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 5)));
            np(6) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 6)));
            
            % Compute the average velocity.
            av_vel = (1/(grain_size.^2))*(np(1)*c1 + np(2)*c2 + np(3)*c3 + np(4)*c4 + np(5)*c5 + np(6)*c6);
            
            % Store the velocity components.
            av_vel_x_comps(currval) = av_vel(1);
            av_vel_y_comps(currval) = av_vel(2);
            
            % Store the positional coordinates.
            av_vel_x_coords(currval) = i;
            av_vel_y_coords(currval) = j;
            
            currval = currval + 1;
        end
    end
    
    quiver(av_vel_x_coords, av_vel_y_coords, av_vel_x_comps, av_vel_y_comps);
    
    hold on;
    plot([1; grain_x], [0.75; 0.75], 'k-');
    hold on;
    plot([1; grain_x], [grain_y + 0.25; grain_y + .25], 'k-');
    
    % Display the flow obstacle.
    obstacle_x = zeros(1, nnz(obstacle));
    obstacle_y = zeros(1, nnz(obstacle));
    k = 1;
    
     for (i = 1:1:numnodes_x)
        for (j = 1:1:numnodes_y)
            if (obstacle(i, j) == 1)
                obstacle_x(k) = 0.5 + (numnodes_x ./ (grain_size .* (numnodes_x - 1))) .* (i - 1);
                obstacle_y(k) = 0.5 + (numnodes_y ./ (grain_size .* (numnodes_y - 1))) .* (j - 1);
                k = k + 1;
            end
        end
    end
    
    hold on; 
    plot(obstacle_x, obstacle_y, 'r-');
    axis equal;
    
end    
