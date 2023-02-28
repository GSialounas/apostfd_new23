classdef dog
    
    properties
        breed 
        colour
        bark_pitch
        size
    end
    
   methods 
       function obj = dog(breed,colour,bark_pitch,size)
           obj.breed = breed;
           obj.colour = colour;
           obj.bark_pitch = bark_pitch;
           obj.size = size;
       end
           
       function identify_dog(obj)
           disp("This is a "+ obj.size + ", " + obj.colour+ " " + obj.breed + ", with a " + obj.bark_pitch+ " bark")
   
       end
   end
end