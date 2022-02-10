function [Reactants_Text, Products_Text] = Coeff2Rxn(xi, xi_prime, Species)
% This function takes the coefficients xi of the reactants, the
% coefficients xi_prime of the products, and the list of species and
% produce the corresponding reaction text.

%% Reactants
xi_binary = xi>0;
if sum(xi_binary) == 0
    Reactants_Text = 'phi';
else
    ReactantSpecies = Species(xi_binary);
    ReactantCoefficients = xi(xi_binary);
    ReactantCoefficients_String = num2str(ReactantCoefficients);
    for i = 1 : length(ReactantCoefficients_String)
        if ReactantCoefficients(i) == 1
            ReactantCoefficients_String(i) = ' ';
        end
    end
    Reactants_Text = [ReactantCoefficients_String(1), ReactantSpecies{1}];
    for i = 2 : length(ReactantSpecies)
        Reactants_Text = [Reactants_Text, ' + ', ReactantCoefficients_String(i), ReactantSpecies{i}];
    end
end

%% Products
xi_binary_prime = xi_prime>0;
if sum(xi_binary_prime) == 0
    Products_Text = 'phi';
else
    ProductSpecies = Species(xi_binary_prime);
    ProductCoefficients = xi_prime(xi_binary_prime);
    ProductCoefficients_String = num2str(ProductCoefficients);
    for i = 1 : length(ProductCoefficients_String)
        if ProductCoefficients(i) == 1
            ProductCoefficients_String(i) = ' ';
        end
    end
    Products_Text = [ProductCoefficients_String(1), ProductSpecies{1}];
    for i = 2 : length(ProductSpecies)
        Products_Text = [Products_Text, ' + ', ProductCoefficients_String(i), ProductSpecies{i}];
    end
end

end

