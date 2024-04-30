function [coeff] = pad_coeff_to_new_min_index(coeff, current_min_index, new_min_index)
    [row, col] = size(coeff);
    num_to_pad = current_min_index - new_min_index;
    top = zeros(num_to_pad, col);
    right = zeros(num_to_pad + row, num_to_pad);
    coeff = [top; coeff];
    coeff = [right, coeff];
end
