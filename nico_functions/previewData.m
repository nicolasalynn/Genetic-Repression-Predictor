function previewData(some_matrix, elements)

[row, col, dim] = size(some_matrix);

    for i = 1:dim
        disp(some_matrix(1:elements, 1:elements, i))
    end

end