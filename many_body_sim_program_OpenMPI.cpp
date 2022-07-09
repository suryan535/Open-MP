/*
						  OPEN MPI Group Assignment
	Group - H
	Members
	Rekha Lokesh       - 19EC10052
	Jothi Prakash      - 19EC30028
	Ashlesh Kumar      - 19EC10080
	Chaitanya Bhargav N - 19EC10016
	Sampara Sai Charan - 19EC10057
*/

#include <iostream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <omp.h>
#include <sstream>
#include <fstream>

using namespace std;

// Cleaning the String to get the required values
string cleanString(string line)
{
	reverse(line.begin(), line.end());
	while (line.back() != '=')
		line.pop_back();
	line.pop_back();
	reverse(line.begin(), line.end());
	return line;
}

int main()
{

	// File stream to read the input and write the output
	fstream fptr;
	fstream res;
	res.open("Coordinates.txt", ios::out);

	fptr.open("Trajectory.txt", ios::in);

	// Specifying the number of threads
	int threadCnt = 6;

	// Handling the file
	if (fptr.is_open())
	{

		// Getting the variable values
		string line;
		double width, length, depth;
		int body_cnt;
		double del_t, mass, radius;

		getline(fptr, line);
		line = cleanString(line);
		width = stod(line);
		getline(fptr, line);
		line = cleanString(line);
		length = stod(line);
		getline(fptr, line);
		line = cleanString(line);
		depth = stod(line);
		getline(fptr, line);
		line = cleanString(line);
		body_cnt = stoi(line);
		getline(fptr, line);
		line = cleanString(line);
		radius = stod(line);
		getline(fptr, line);
		line = cleanString(line);
		mass = stod(line);
		getline(fptr, line);
		line = cleanString(line);
		del_t = stod(line);

		getline(fptr, line);

		double coord[body_cnt][3];
		int i = 0;
		while (getline(fptr, line))
		{
			stringstream ss(line);
			ss >> coord[i][0] >> coord[i][1] >> coord[i][2];
			i++;
		}

		// Closing the file after reading
		fptr.close();

		// Other needed values
		double del_t_by_m = del_t / mass;
		double minDis = 2 * radius;
		int totalTime = 1000;
		double force[body_cnt][3];
		double velocity[body_cnt][3];
		memset(velocity, 0, sizeof(velocity));

		// Storing the values in the output file
		res << "Width: " << width << "\n";
		res << "Length: " << length << "\n";
		res << "Depth: " << depth << "\n";
		res << "Body Count: " << body_cnt << "\n";
		res << "Body radius: " << radius << "\n";
		res << "Time Step: " << del_t << "\n";
		res << "Iterations: " << 7200 << "\n";

		// Measuring the time
		double start;
		double end;
		start = omp_get_wtime();

		// Code for the simulation
		for (int n = 0; n < totalTime; n++)
		{

			// Finding the force on each body due to interbody forces
#pragma omp parallel for num_threads(threadCnt) collapse(2)
			for (int i = 0; i < body_cnt; i++)
			{
				for (int j = 0; j < body_cnt; j++)
				{
					if (i == j)
						continue;
					double dis = coord[j][0] - coord[i][0];
					if (dis >= 0.0)
					{
						dis = max(dis, minDis);
					}
					else
					{
						dis = min(dis, -minDis);
					}
					force[i][0] = (mass * mass) / (dis * dis);
					dis = coord[j][1] - coord[i][1];
					if (dis >= 0.0)
					{
						dis = max(dis, minDis);
					}
					else
					{
						dis = min(dis, -minDis);
					}
					force[i][1] = (mass * mass) / (dis * dis);
					dis = coord[j][2] - coord[i][2];
					if (dis >= 0.0)
					{
						dis = max(dis, minDis);
					}
					else
					{
						dis = min(dis, -minDis);
					}
					force[i][2] = (mass * mass) / (dis * dis);
				}
			}

			// Changing the velocity and coordinates based on the forces
#pragma omp parallel for num_threads(threadCnt)
			for (int i = 0; i < body_cnt; i++)
			{

				velocity[i][0] += del_t_by_m * force[i][0] * 0.5;
				velocity[i][1] += del_t_by_m * force[i][1] * 0.5;
				velocity[i][2] += del_t_by_m * force[i][2] * 0.5;

				coord[i][0] += velocity[i][0] * del_t;
				coord[i][1] += velocity[i][1] * del_t;
				coord[i][2] += velocity[i][2] * del_t;

				velocity[i][0] += del_t_by_m * force[i][0] * 0.5;
				velocity[i][1] += del_t_by_m * force[i][1] * 0.5;
				velocity[i][2] += del_t_by_m * force[i][2] * 0.5;

				if (coord[i][0] > width || coord[i][0] < 0)
				{
					velocity[i][0] = -velocity[i][0];
					coord[i][0] = min(coord[i][0], width);
					coord[i][0] = max(coord[i][0], 0.0);
				}

				if (coord[i][1] > length || coord[i][1] < 0)
				{
					velocity[i][1] = -velocity[i][1];
					coord[i][1] = min(coord[i][1], length);
					coord[i][1] = max(coord[i][1], 0.0);
				}

				if (coord[i][2] > depth || coord[i][2] < 0)
				{
					velocity[i][2] = -velocity[i][2];
					coord[i][2] = min(coord[i][2], depth);
					coord[i][2] = max(coord[i][2], 0.0);
				}
			}

			// Checking for collisions and applying Newton physics for momentum
#pragma omp parallel for num_threads(threadCnt) collapse(2)
			for (int i = 0; i < body_cnt; i++)
			{
				for (int j = 0; j < body_cnt; j++)
				{
					if (i == j)
						continue;

					double disX = abs(coord[i][0] - coord[j][0]);
					disX *= disX;
					double disY = abs(coord[i][1] - coord[j][1]);
					disY *= disY;
					double disZ = abs(coord[j][2] - coord[i][2]);
					disZ *= disZ;

					double dis = disX + disY + disZ;
					dis = sqrt(dis);
					if (dis < minDis && i < j)
					{

						swap(velocity[i][0], velocity[j][0]);

						swap(velocity[i][1], velocity[j][1]);

						swap(velocity[i][2], velocity[j][2]);
					}
				}
			}

			// Displaying progress after every 100 iterations
			if (n % 100 == 0)
			{
				cout << n << "\n";
				res << "Iteration : " << n / 100 << "\n";
				for (int j = 0; j < body_cnt; j++)
				{
					res << coord[j][0] << " " << coord[j][1] << " " << coord[j][2] << "\n";
				}
			}
		}

		// Closing the output file
		res.close();

		// Displaying the total elapsed time
		end = omp_get_wtime();
		printf("Work took %f seconds\n", end - start);
	}
	else
	{
		// In case of any error
		cout << "Error in Opening the File\n";
	}

	return 0;
}