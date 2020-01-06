function reconstructed_data = reconstruct_data(deconstructed, blueprint)


    reconstructed_data = zeros(1, 5805);
    
    for i = 1:length(blueprint)

        reconstructed_data(1, blueprint(i)) = deconstructed(i);

    end


end