function indices = block_wave_edges(block_wave)
    % This function finds the indices indicating the start and end of each
    % block in a block wave.

    % Input arguments
    % block_wave: you know, a block wave, between 0 and 1

    % Output
    % Indices: pairs of indices for each block, (n_blocks, 2)

    indices = [];
    high = false;

    for i=1:size(block_wave, 1)
        if not(high) && block_wave(i)
            high = true;
            index = i;
        end
        if high && not(block_wave(i))
            high = false;
            indices = [indices; [index, i-1]];
        end
        if high && (i == size(block_wave, 1))
            indices = [indices; [index, i]];
        end
    end
end
