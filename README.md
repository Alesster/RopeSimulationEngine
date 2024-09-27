Here We consider the rope simulation. It is based on position based dynamics. The main reference for simulation is taken from:
https://www.owlree.blog/posts/simulating-a-rope.html
The main idea: 
 - The rope is devided into N elements. 
 - Each element is simulated by using the Equations of Motion:
 - Euler method or Verlet method.
 - Additionally the Jakobsen Method (Enforcing Constraints) is used for the rope relaxation
