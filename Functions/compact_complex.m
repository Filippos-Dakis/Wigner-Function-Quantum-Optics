function str = compact_complex(z)
    real_part = real(z);
    imag_part = imag(z);
    
    % Format real part
    if abs(real_part - round(real_part)) < eps
        real_str = num2str(round(real_part));
    else
        real_str = num2str(real_part);
    end
    
    % Format imaginary part
    if abs(imag_part - round(imag_part)) < eps
        imag_str = [num2str(round(imag_part)) 'i'];
    else
        imag_str = [num2str(imag_part) 'i'];
    end

    
    % Combine real and imaginary parts
    if real_part == 0 && imag_part == 0
        str = '0';
    elseif abs(imag_part) <1e-15
        str = real_str;
    elseif abs(real_part) < 1e-15
        str = imag_str;
    elseif imag_part > 0
        str = [real_str ' + ' imag_str];
    else
        str = [real_str ' - ' num2str(-imag_part) 'i'];
    end
end
