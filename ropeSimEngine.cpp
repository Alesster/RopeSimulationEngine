#include<iostream>
#include<math.h>
#include<assert.h>
#include <iomanip>
#include <string>
#include<fstream>

using namespace std;
typedef float f;

float const PI = 3.14159;

// ========== BEGIN CLASS vector3d ==========

class vector3d
{
public:
    f x, y, z;
    vector3d()  //constructor
    {
        x = 0;
        y = 0;
        z = 0;
    }

    vector3d(f x1, f y1, f z1 = 0)     //initializing object with values.
    {
        x = x1;
        y = y1;
        z = z1;
    }
    vector3d(const vector3d& vec);    //copy constructor
    vector3d operator+(const vector3d& vec);    //addition
    vector3d& operator+=(const vector3d& vec);  ////assigning new result to the vector
    vector3d operator-(const vector3d& vec);    //substraction
    vector3d operator-();
    vector3d& operator-=(const vector3d& vec);  //assigning new result to the vector
    vector3d operator*(f value);    //multiplication
    vector3d& operator*=(f value);  //assigning new result to the vector.
    vector3d operator/(f value);    //division
    vector3d& operator/=(f value);  //assigning new result to the vector
    vector3d& operator=(const vector3d& vec);
    f dot_product(const vector3d& vec); //scalar dot_product
    vector3d cross_product(const vector3d& vec);    //cross_product
    f magnitude();  //magnitude of the vector
    vector3d normalization();   //nor,malized vector
    f square(); //gives square of the vector

    f distance(const vector3d& vec);    //gives distance between two vectors
    f show_X(); //return x
    f show_Y(); //return y
    f show_Z(); //return z
    void disp();    //display value of vectors
};

vector3d::vector3d(const vector3d& vec)
{
    x = vec.x;
    y = vec.y;
    z = vec.z;
}

//addition

vector3d vector3d ::operator+(const vector3d& vec)
{
    return vector3d(x + vec.x, y + vec.y, z + vec.z);
}

vector3d& vector3d ::operator+=(const vector3d& vec)
{
    x += vec.x;
    y += vec.y;
    z += vec.z;
    return *this;
}
//substraction//
vector3d vector3d ::operator-(const vector3d& vec)
{
    return vector3d(x - vec.x, y - vec.y, z - vec.z);
    //return vector3d(-vec.x, -vec.y, -vec.z);
}

vector3d vector3d ::operator-()
{
    //return vector3d(x - vec.x, y - vec.y, z - vec.z);
    return vector3d(-x, -y, -z);
}

vector3d& vector3d::operator-=(const vector3d& vec)
{
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
    return *this;
}

//scalar multiplication

vector3d vector3d ::operator*(f value)
{
    return vector3d(x * value, y * value, z * value);
}



vector3d& vector3d::operator*=(f value)
{
    x *= value;
    y *= value;
    z *= value;
    return *this;
}

//scalar division
vector3d vector3d ::operator/(f value)
{
    assert(value != 0);
    return vector3d(x / value, y / value, z / value);
}

vector3d& vector3d ::operator/=(f value)
{
    assert(value != 0);
    x /= value;
    y /= value;
    z /= value;
    return *this;
}


vector3d& vector3d::operator=(const vector3d& vec)
{
    x = vec.x;
    y = vec.y;
    z = vec.z;
    return *this;
}

//Dot product
f vector3d::dot_product(const vector3d& vec)
{
    return x * vec.x + vec.y * y + vec.z * z;
}

//cross product
vector3d vector3d::cross_product(const vector3d& vec)
{
    f ni = y * vec.z - z * vec.y;
    f nj = z * vec.x - x * vec.z;
    f nk = x * vec.y - y * vec.x;
    return vector3d(ni, nj, nk);
}

f vector3d::magnitude()
{
    return sqrt(square());
}

f vector3d::square()
{
    return x * x + y * y + z * z;
}

vector3d vector3d::normalization()
{
    assert(magnitude() != 0);
    *this /= magnitude();
    return *this;
}

f vector3d::distance(const vector3d& vec)
{
    vector3d dist = *this - vec;
    return dist.magnitude();
}

f vector3d::show_X()
{
    return x;
}

f vector3d::show_Y()
{
    return y;
}

f vector3d::show_Z()
{
    return z;
}

void vector3d::disp()
{
    cout << x << " " << y << " " << z << endl;
}

//  ========== END OF CLASS vector3d ==========

//  ========== BEGIN CLASS Mass ==========
class Mass
{
public:

    float     m;        // The value of mass
    vector3d  pos;      // Current Position
    vector3d  oldpos;   // Old position
    vector3d  temp;     // Temporal value
    vector3d  vel;      // Velocity
    vector3d  force;    // Acting force

    Mass(float m)       // Constructor
    {
        this->m = m;
    }

    void applyForce(vector3d force)
    {
        this->force += force;    // Adding external force
    }

    void init()          // Force deleting
    {
        force.x = 0;
        force.y = 0;
        force.z = 0;
    }

    void simulate(float dt)
    {
        vel += (force / m) * dt;
        //pos += vel * dt;
        temp = pos;
        pos += pos - oldpos + (force / m) * dt * dt;
        oldpos = temp;
    }
};
// ========== END OF CLASS Mass ==========

//  ========== BEGIN CLASS Simulation ==========

class Simulation
{
public:
    int numOfMasses;
    Mass** masses;

    Simulation(int numOfMasses, float m)
    {
        this->numOfMasses = numOfMasses;
        masses = new Mass * [numOfMasses];

        for (int i = 0; i < numOfMasses; i++)
        {
            masses[i] = new Mass(m);
        }
    }

    virtual void release()
    {
        for (int i = 0; i < numOfMasses; i++)
        {
            delete[] masses[i];
            masses[i] = NULL;
        }

        delete[] masses;
        masses = NULL;
    }

    Mass* getMass(int index)
    {
        if (index < 0 || index >= numOfMasses) return NULL;

        return masses[index];
    }

    virtual void init()
    {
        for (int i = 0; i < numOfMasses; i++)
        {
            masses[i]->init();
        }
    }

    virtual void solve()
    {

    }

    virtual void simulate(float dt)
    {
        for (int i = 0; i < numOfMasses; i++)
        {
            masses[i]->simulate(dt);
        }
    }

    virtual void operate(float dt)
    {
        init();
        solve();
        simulate(dt);
    }
};
//  ========== END OF CLASS Simulation ==========

// ========== class SPRING ==========

class Spring
{
public:
    Mass* mass1;
    Mass* mass2;
    vector3d gravitation;

    float springConstant;       // Spring constant
    float springLength;         // The length without any force
    float frictionConstant;     // Friction constant

    Spring(Mass* mass1, Mass* mass2, float springConstant,
        float springLength, float frictionConstant)
    {
        this->springConstant = springConstant;
        this->springLength = springLength;
        this->frictionConstant = frictionConstant;

        this->mass1 = mass1;
        this->mass2 = mass2;

        this->gravitation = gravitation;
    }

    void solve()
    {
        vector3d springVector = mass1->pos - mass2->pos;
        float r = springVector.magnitude();

        vector3d force;
        if (r != 0)
        {
            force += (springVector / r) * (r - springLength) * (-springConstant);
        }
        force += -(springVector / r) * (springVector / r).dot_product(mass1->vel - mass2->vel) * frictionConstant;
        mass1->applyForce(force);
        mass2->applyForce(-force);
    }
};

// ========== END class SPRING ==========

// ========== CLASS RopeSimulation  ==========

class RopeSimulation : public Simulation
{
public:
    Spring** springs; // Springs connect numOfMasses of masses.
    vector3d gravitation; // Gravitational acceleration (Applied to all masses)
    vector3d ropeConnectionPos; // First point of Rope (with index 0) that used to set the position
    vector3d ropeConnectionVel; // Velocity of moving of ropeConnectionPos. It is possible to swing the Rope
    float groundRepulsionConstant; // Repultion constant of ground. How much the ground is repulsing the Rope.
    float groundFrictionConstant; // Friction constant of the ground for sliding the Rope on the ground.    
    float groundAbsorptionConstant; // Absorption constant of the ground. Used for the vertical collisions of the Rope with the ground.
    float groundHeight; // Y coordinate of the ground.
    float airFrictionConstant; // Friction constant of the air. 

    // Constructor with 11 parameters:
    RopeSimulation(                  // Constructor
        int numOfMasses,               // 1. Number of masses
        float m,                       // 2. Weight of each mass
        float springConstant,          // 3. How tight is the spring
        float springLength,            // 4. Length of relaxed spring
        float springFrictionConstant,  // 5. Internal friction
        vector3d gravitation,          // 6. Gravtation acceleration
        float airFrictionConstant,     // 7. Air friction
        float groundRepulsionConstant, // 8. Ground repulsion
        float groundFrictionConstant,  // 9. Ground friction
        float groundAbsorptionConstant,// 10. Ground absorption
        float groundHeight             // 11. Height of the ground
    ) : Simulation(numOfMasses, m) // The class that create masses with weight m. 
    {
        this->gravitation = gravitation;
        this->airFrictionConstant = airFrictionConstant;
        this->groundFrictionConstant = groundFrictionConstant;
        this->groundRepulsionConstant = groundRepulsionConstant;
        this->groundAbsorptionConstant = groundAbsorptionConstant;
        this->groundHeight = groundHeight;

        for (int i = 0; i < numOfMasses; i++)
        {// INITIAL POSITION:
            masses[i]->pos.x = i * springLength; // *sin(PI / 4);
            masses[i]->pos.y = 0.0; // -i * springLength; // *sin(PI / 4);
            masses[i]->pos.z = 0;

            masses[i]->oldpos = masses[i]->pos;
        }

        springs = new Spring * [numOfMasses - 1];

        for (int i = 0; i < numOfMasses - 1; i++)
        {
            springs[i] = new Spring(masses[i], masses[i + 1],
                springConstant, springLength, springFrictionConstant);
        }
    }

    // Adding all forces
    void solve() override
    {
        for (int i = 0; i < numOfMasses - 1; i++)
        {
            springs[i]->solve();
        }

        for (int i = 0; i < numOfMasses; i++)
        {
            masses[i]->applyForce(gravitation * masses[i]->m);
            masses[i]->applyForce(-masses[i]->vel * airFrictionConstant);

            if (masses[i]->pos.y < groundHeight)
            {
                vector3d v;

                v = masses[i]->vel;
                v.y = 0;

                // Friction force of the ground 
                masses[i]->applyForce(-v * groundFrictionConstant);

                v = masses[i]->vel;
                v.x = 0;
                v.z = 0;

                // Absorption and repulsion by the ground
                if (v.y < 0)
                {
                    masses[i]->applyForce(-v * groundAbsorptionConstant);

                    vector3d force = vector3d(0, groundRepulsionConstant, 0) *
                        (groundHeight - masses[i]->pos.y);
                    masses[i]->applyForce(force);
                }
            }
        }
    }

    // Rope first point movement simulation
    void simulate(float dt) override
    {
        Simulation::simulate(dt);

        ropeConnectionPos += ropeConnectionVel * dt;

        if (ropeConnectionPos.y < groundHeight)
        {
            ropeConnectionPos.y = groundHeight;
            ropeConnectionVel.y = 0;
        }

        masses[0]->pos = ropeConnectionPos;
        masses[0]->vel = ropeConnectionVel;
    }

    // Rope first point velocity
    void setRopeConnectionVel(vector3d ropeConnectionVel)
    {
        this->ropeConnectionVel = ropeConnectionVel;
    }
};

// ========== END CLASS RopeSimulation  ==========

int main()
{   
    int frame = 0;          // Temporary value of the frame to show position
    int numOfMasses = 21;   // Number of masses
    float dt = 0.001f;      // Time interval
    float g = 9.81f;        // Gravitational acceleration
    float mass = 0.05f;     // Mass of each particle
    float length = 0.05f;   // Length of relaxed spring of each link
    float delta = 0.001;    // 1.0 mm - How much the first link is streching under the gravitation.
    float springConstant = (numOfMasses - 1) * mass * g / delta; // Spring Constant
    float springFrictionConstant = 2.0;     // Constant of internal friction
    float airFrictionConstant = 0.02;       // Air friction
    float groundRepulsionConstant = 100.0f; // Ground repulsion
    float groundFrictionConstant = 0.2f;    // Ground friction
    float groundAbsorptionConstant = 2.0f;  // Ground absorption
    float groundHeight = -1.5f;             // Height of the ground

    vector3d ropeConnectionVel = vector3d(0, 0, 0);  // Rope first point velocity

    RopeSimulation* ropeSimulation =
        new RopeSimulation(
            numOfMasses,                // Number of masses
            mass,                       // Mass of each particle
            springConstant,             // Spring Constant
            length,                     // Length of relaxed spring
            springFrictionConstant,     // Constant of internal friction
            vector3d(0, -g, 0),         // Gravitational acceleration
            airFrictionConstant,        // Air friction constant
            groundRepulsionConstant,    // Constant of ground repulsion
            groundFrictionConstant,     // Constant of ground friction
            groundAbsorptionConstant,   // Constant of ground absorption
            groundHeight);              // Height of ground

    ropeSimulation->setRopeConnectionVel(ropeConnectionVel);
    for (int i = 0; i < 1000; i++)
    {
        /*if (i == 100) {
            ropeConnectionVel = vector3d(0.0, 0, 0);
            ropeSimulation->setRopeConnectionVel(ropeConnectionVel);
        }*/
       
        ropeSimulation->operate(dt);

        if (i == frame) {
            ofstream fout(string("output_" + to_string(i) + ".dat").c_str());
            for (int j = 0; j < numOfMasses; j++)
            {
                vector3d pos = ropeSimulation->masses[j]->pos;  
                vector3d vel = ropeSimulation->masses[j]->vel;              
                //pos.disp();
                fout << fixed << setprecision(6) << pos.x << "\t " << pos.y << endl;
            }
            fout.close();
            frame += 30;
        }
    }

    ropeSimulation->release();

    return 0;
}