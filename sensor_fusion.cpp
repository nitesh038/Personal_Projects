/* Copyright (2018) KPIT Technologies */
#include <cstdio>
#include"math.h"
#include"base\ReadCamParams.h"
#include"base\ReadCSV_VED_CSV.h"
#include"opencv2/opencv.hpp"
#include"opencv2/core/core.hpp"
#define NUMBER_OF_ATTRIBUTES_IN_PED 39
#define NUMBER_OF_ATTRIBUTES_IN_VED 25
int flag_for_PED = -1;
int flag_for_VED = -1;
int flag_for_CSV = 1;
std::string Main_Csv;
std::string FISHEYE_Csv;
std::string TELE_Csv;
std::string ConfigFile;
int image_width = 1280;
int image_height = 964;
int image_height_TELE = 1084;
enum Camera {
  STANDARD = 0,
  FISHEYE = 1,
  TELE = 2
}Cam;
VED_Header H_main;
VED_Header H_fish_eye;
VED_Header H_TELE;
struct Points {
  double left;
  double bottom;
  double right;
  double top;
  double xpos;
  double ybottom;
  double ytop;
};
Points P;
double f_fe;
double f_std;
double f_TELE;
double no_of_rows_in_std;
double no_of_rows_in_TELE;
double no_of_rows_in_fish_eye;
int num_of_rows_main = 0;
int num_of_rows_fish_eye = 0;
int num_of_rows_TELE = 0;
double Global_Start_Ind;
int Global_Image_St_Row;
std::vector<double>grab_index_main;
std::vector<double>grab_index_sorted_main;
std::vector<double>grab_index_fish_eye;
std::vector<double>grab_index_TELE;
std::vector<double>all_distance;
std::vector<double>grab_index_only_in_a;
std::vector<double>row_index_only_in_FISHEYE;
std::vector<double>row_index_only_in_TELE;
std::map<double, double>fish_eye_row_mapped_to_grab_index;
std::map<double, double>TELE_row_mapped_to_grab_index;
std::map<double, double>main_row_mapped_to_grab_index;
std::map<double, double>distance_mapped_to_grab_index;
std::map<double, double>fish_eye_grab_index_mapped_to_row;
std::map<double, double>TELE_grab_index_mapped_to_row;
std::map<double, double>main_grab_index_mapped_to_row;
std::map<double, Points>fish_eye_coordinates_mapped_to_fusion_row;
std::map<double, Points>TELE_coordinates_mapped_to_fusion_row;
std::map<double, Points>main_coordinates_mapped_to_grab_index;
std::map<double, Points>sorted_main_coordinates_mapped_to_grab_index;
std::map<double, Points>fish_eye_coordinates_all_rows;
std::map<double, Points>TELE_coordinates_all_rows;
// final mapping to consider to transform
std::map<double, Points>transformed_fish_eye_coordinates_mapped_to_fusion_row;
std::map<double, Points>transformed_TELE_coordinates_mapped_to_fusion_row;
std::vector<Points>Fish_Eye_Points_in_each_row_for_fusion;
std::vector<Points>TELE_Points_in_each_row_for_fusion;
std::vector<double>X_value_only_in_FISHEYE_to_be_modified;
std::vector<double>Y_value_only_in_FISHEYE_to_be_modified;
// this will map row index of F.E to distance and also helps
// in discarding row index whose grab index not present in distance file
std::map<double, double>fish_eye_row_index_mapped_to_distance;
// this will map row index of F.E to distance and also helps
// in discarding row index whose grab index not present in distance file
std::map<double, double>TELE_row_index_mapped_to_distance;
// Hold Data for main,FISHEYE, TELE output
std::map<double, VED_Header>map_main_output_to_row;
std::map<double, VED_Header>map_FISHEYE_output_to_row;
std::map<double, VED_Header>map_TELE_output_to_row;
// for final fused output
std::map<double, VED_Header>map_main_output_to_GrabIndex;
std::map<double, VED_Header>map_FISHEYE_output_to_GrabIndex;
std::map<double, VED_Header>map_TELE_output_to_GrabIndex;
std::map<double, VED_Header>map_final_main_output_to_GrabIndex;
std::map<double, VED_Header>map_final_main_output_to_row_number;
std::map<double, VED_Header>map_SortOutput_output_to_GrabIndex;
std::map<double, double>Corrected_distance_mapped_to_grab_index;
std::map<double, double>Sorted_distance_mapped_to_grab_index;
std::map<double, double>::iterator Dist_Iterator;
std::map<double, double>main_xpos_mapped_to_grab_index;
std::map<double, double>fish_eye_xpos_mapped_to_grab_index;
std::map<double, double>TELE_xpos_mapped_to_grab_index;
std::map<int, int>main_object_type_mapped_to_grab_index;
std::map<int, int>fish_eye_object_type_mapped_to_grab_index;
std::map<int, int>TELE_object_type_mapped_to_grab_index;
std::map<int, int>::iterator main_object_type_iterator;
std::map<int, int>::iterator fish_eye_object_type_iterator;
std::map<int, int>::iterator TELE_object_type_iterator;
std::vector<VED_Header>H_to_write;
struct image_coordinates_in_X_Y {
  int X;
  int Y;
};
struct point_in_image_coordinates {
  image_coordinates_in_X_Y topLeft;
  image_coordinates_in_X_Y topRight;
  image_coordinates_in_X_Y bottomLeft;
  image_coordinates_in_X_Y bottomRight;
  image_coordinates_in_X_Y  xpos_ybottom;
  image_coordinates_in_X_Y  xpos_ytop;
};
point_in_image_coordinates points_;
double radial_distance_main(double r_d);
double calculate_r_d(int type_of_Camera, int X, int Y);
void new_undistorted_points(int type_of_Camera, double r_u, double r_d, 
  double X_d, double Y_d, int *X_u, int *Y_u);
std::vector<double> Gl_SortedGrab;
// Read the output file and rewrite it based on sorted GrabIndex
void sorted_output() {
  int i = 0;
  int num_main_rows = 0;
  int num_FISHEYE_rows = 0;
  int num_TELE_rows = 0;
  int num_rows = 0;
  int Row_Start = 0;
  double value_to_insert;
  std::map<double, VED_Header>::iterator itr_Cam;
  double min_grab_id = 99999999;
  double max_grab_id = -1;
  std::ifstream file("Output.csv");
  std::ofstream outfile;
  outfile.open("SortedOutput.csv", std::ofstream::out|std::ofstream::trunc);
  CSVRow row;
  while(file >> row) {
  H_main.Keyframe = row[0];
	H_main.VCL = row[1];
	H_main.Perfect = row[2];
	H_main.camPort = row[3];
	H_main.originalID	= row[4];
	H_main.left = row[5];
	H_main.bottom	=	row[6];
	H_main.right = row[7];
	H_main.top = row[8];
	H_main.xpos	=	row[9];
	H_main.ybottom = row[10];
	H_main.ytop	=	row[11];
	H_main.type	=	row[12];
	H_main.truckType = row[13];
	H_main.lane =	row[14];
	H_main.oncoming	=	row[15];
	H_main.CIPV	=	row[16];
	H_main.mobility	=	row[17];
	H_main.occluded	=	row[18];
	H_main.blinking	=	row[19];
	H_main.braking	=	row[20];
	H_main.grabIndex = row[21];
	H_main.enterHostLane = row[22];
	H_main.rearOcclusion =	row[23];
	H_main.wheelVisibility = row[24];
	double grabIndex = atoi(H_main.grabIndex.c_str());

		if(grabIndex == 0)
			continue;

		num_rows++;
		double camPort = atoi(H_main.camPort.c_str());
		// Set CamPort for Fusion to 13
		H_main.camPort = "13";

		if(grabIndex < min_grab_id)
			min_grab_id = grabIndex;
		if(grabIndex > max_grab_id)
			max_grab_id = grabIndex;

		if(camPort == 0)
		{	// Main
			map_main_output_to_row.insert(std::pair<double,VED_Header> (grabIndex, H_main));
			num_main_rows++;
		}else if(camPort == 2)
		{	//Fish Eye
			map_FISHEYE_output_to_row.insert(std::pair<double,VED_Header> (grabIndex, H_main));
			num_FISHEYE_rows++;
		}else if(camPort == 1)
		{	// TELE
			map_TELE_output_to_row.insert(std::pair<double,VED_Header> (grabIndex, H_main));
			num_TELE_rows++;
		}
		else
		{	
			std::cout<<"Invalid Cam port"<<std::endl;
		}
    }

	if(num_rows > 0)
	{
		for(int grab_id = min_grab_id; grab_id <= max_grab_id; grab_id++)
		{

			if(map_main_output_to_row.count(grab_id) > 0)
			{
				itr_Cam = map_main_output_to_row.find(grab_id);

			}else if(map_TELE_output_to_row.count(grab_id) > 0)
			{
				itr_Cam = map_TELE_output_to_row.find(grab_id);

			}else if(map_FISHEYE_output_to_row.count(grab_id) > 0)
			{
				itr_Cam = map_FISHEYE_output_to_row.find(grab_id);
			}
			else
			{	continue;
			}

			H_main = itr_Cam->second;
			char Conv_Row[10] = {'\0'};
			sprintf(Conv_Row,"%d",Row_Start);
			std::string Str1(Conv_Row);
			Row_Start++;
			value_to_insert = atoi(H_main.grabIndex.c_str()) ;
			Gl_SortedGrab.push_back(value_to_insert);
			// H_main.Keyframe			 	=				Str1;
			map_SortOutput_output_to_GrabIndex.insert(std::pair<double, VED_Header>
			(value_to_insert,H_main));

			// outfile.write
      outfile << H_main.Keyframe << ",";		
			outfile << H_main.VCL	<< ",";		
			outfile << H_main.Perfect	<< ",";		
			outfile << H_main.camPort	<< ",";		
			outfile << H_main.originalID << ",";	
			outfile << H_main.left << ",";	
			outfile << H_main.bottom << ",";	
			outfile << H_main.right << ",";	
			outfile << H_main.top << ",";	
			outfile << H_main.xpos << ",";	
			outfile << H_main.ybottom	<< ",";	
			outfile << H_main.ytop << ",";	
			outfile << H_main.type << ",";	
			outfile << H_main.truckType	<< ",";	
			outfile << H_main.lane << ",";	
			outfile << H_main.oncoming	<< ",";	
			outfile << H_main.CIPV << ",";	
			outfile << H_main.mobility << ",";	
			outfile << H_main.occluded << ",";	
			outfile << H_main.blinking	<< ",";	
			outfile << H_main.braking	<< ",";	
			outfile << H_main.grabIndex << ",";	
			outfile << H_main.enterHostLane	<< ",";	
			outfile << H_main.rearOcclusion << ",";	
			outfile << H_main.wheelVisibility;	
			outfile << std::endl;					

		}// for

	}// if

	outfile.close();
}

void Fused_Write();
void WriteToFused(VED_Header* He_Data);

void Fused_Write()
{
	int Total_Index = Gl_SortedGrab.size();
	int Ind;
	VED_Header  H_For_Fus;
  for(Ind = 0; Ind < Total_Index ; Ind++)
	{
		if(map_main_output_to_GrabIndex.count(Gl_SortedGrab[Ind]) > 0)
		{
			H_For_Fus = map_main_output_to_GrabIndex[Gl_SortedGrab[Ind]];
			WriteToFused(&H_For_Fus);
		}
		if(map_TELE_output_to_GrabIndex.count(Gl_SortedGrab[Ind]) > 0)
		{
			H_For_Fus = map_TELE_output_to_GrabIndex[Gl_SortedGrab[Ind]];
			WriteToFused(&H_For_Fus);
		}
		if(map_FISHEYE_output_to_GrabIndex.count(Gl_SortedGrab[Ind]) > 0)
		{
			H_For_Fus = map_FISHEYE_output_to_GrabIndex[Gl_SortedGrab[Ind]];
			WriteToFused(&H_For_Fus);
		}
		if(map_SortOutput_output_to_GrabIndex.count(Gl_SortedGrab[Ind]) > 0)
		{
			H_For_Fus = map_SortOutput_output_to_GrabIndex[Gl_SortedGrab[Ind]];
			WriteToFused(&H_For_Fus);
		}

	}

	
}

void WriteToFused(VED_Header* He_Data)
{
		std::ofstream outfile;
		outfile.open ("AllFused_Output.csv" , std::ofstream::out | std::ofstream::app );
    	outfile <<	He_Data->Keyframe	<< ",";		
			outfile <<	He_Data->VCL << ",";		
			outfile <<	He_Data->Perfect << ",";		
			outfile <<	He_Data->camPort << ",";		
			outfile <<	He_Data->originalID	<< ",";	
			outfile <<	He_Data->left	<< ",";	
			outfile <<	He_Data->bottom	<< ",";	
			outfile <<	He_Data->right << ",";	
			outfile <<	He_Data->top << ",";	
			outfile <<	He_Data->xpos	<< ",";	
			outfile <<	He_Data->ybottom << ",";	
			outfile <<	He_Data->ytop	<< ",";	
			outfile <<	He_Data->type	<< ",";	
			outfile <<	He_Data->truckType << ",";	
			outfile <<	He_Data->lane	<< ",";	
			outfile <<	He_Data->oncoming	<< ",";	
			outfile <<	He_Data->CIPV	<< ",";	
			outfile <<	He_Data->mobility	<< ",";	
			outfile <<	He_Data->occluded	<< ",";	
			outfile <<	He_Data->blinking	<< ",";	
			outfile <<	He_Data->braking	<< ",";	
			outfile <<	He_Data->grabIndex << ",";	
			outfile <<	He_Data->enterHostLane << ",";	
			outfile <<	He_Data->rearOcclusion	<< ",";	
			outfile <<	He_Data->wheelVisibility;	
			outfile << std::endl;					
      outfile.close();

}




void cartesian_to_image(int type_of_Camera, image_coordinates_in_X_Y  &X1 , double cart_X , double cart_Y)
{
	int image_width_Cam = image_width;
	int image_height_Cam;
	if(type_of_Camera == TELE )
		 image_height_Cam = image_height_TELE;
	else
		 image_height_Cam = image_height;

	cart_X *= 4;
	cart_Y *= 4;

	cart_X  = cart_X + stCamParams.CamParams[type_of_Camera].Yaw;
	cart_Y =  cart_Y + stCamParams.CamParams[type_of_Camera].Horizon;

	X1.X =  cart_X + image_width_Cam/2;
	X1.Y = -cart_Y + image_height_Cam/2;
}

void image_coordinates_to_cartesian(int type_of_Camera, image_coordinates_in_X_Y &X1 , double img_X , double img_Y)
{
	int image_width_Cam = image_width;
	int image_height_Cam;
	if(type_of_Camera == TELE )
		 image_height_Cam = image_height_TELE;
	else
		 image_height_Cam = image_height;

	img_X =  img_X - image_width_Cam/2;
	img_Y = -img_Y + image_height_Cam/2;

	img_X  = img_X - stCamParams.CamParams[type_of_Camera].Yaw;
	img_Y =  img_Y - stCamParams.CamParams[type_of_Camera].Horizon;

	X1.X = img_X / 4;
	X1.Y = img_Y / 4;
}

void write_to_particular_row(Camera type_of_Camera)
{
	if(type_of_Camera == FISHEYE)
	{
		std::map<double, Points>::iterator itr;
		std::ifstream       file(FISHEYE_Csv);
		CSVRow              row;
		int i = 0 ;
		int no_of_rows_reached = 0;
		int Fish_Eye_Points_in_each_row_counter = 0;
		while(file >> row)
		{	
			if(row_index_only_in_FISHEYE.size()==0)
			{
				break;
			}

			if(i >= row_index_only_in_FISHEYE.size())
			{	break;
			}

			if(no_of_rows_reached == row_index_only_in_FISHEYE.at(i))
			{
		
				if(transformed_fish_eye_coordinates_mapped_to_fusion_row.count(row_index_only_in_FISHEYE[i]) == 0)
				{	
					i++;
					no_of_rows_reached++;
					continue;
				}

				itr = transformed_fish_eye_coordinates_mapped_to_fusion_row.find(row_index_only_in_FISHEYE[i++]);
				std::ofstream myfile;
				myfile.open ("Output.csv" , std::ofstream::out | std::ofstream::app );
				if (!myfile.is_open())
				{
					std::cout << "Output operation not successfully performed\n";
					myfile.close();
					//return -1;
				}

				myfile << row[0] << ",";		
				myfile << row[1] << ",";		
				myfile << row[2] << ",";		
				myfile << row[3] << ",";		
				myfile << row[4] << ",";	
				myfile << (itr->second).left << ",";	
				myfile << (itr->second).bottom << ",";	
				myfile << (itr->second).right << ",";	
				myfile << (itr->second).top << ",";	
				myfile << (itr->second).xpos << ",";	
				myfile << (itr->second).ybottom << ",";	
				myfile << (itr->second).ytop << ",";	
				myfile << row[12] << ",";	
				myfile << row[13] << ",";	
				myfile << row[14] << ",";	
				myfile << row[15] << ",";	
				myfile << row[16] << ",";	
				myfile << row[17] << ",";	
				myfile << row[18] << ",";	
				myfile << row[19] << ",";	
				myfile << row[20] << ",";	
				myfile << row[21] << ",";	// grab-index
				myfile << row[22] << ",";	
				myfile << row[23] << ",";	
				myfile << row[24] << ",";	
				myfile << std::endl;					
				myfile.close();	
						
			}
			no_of_rows_reached++;
		}
	}
	if(type_of_Camera == TELE)
	{
		std::map<double, Points>::iterator itr;
		std::ifstream file(TELE_Csv);
		CSVRow row;
		int i = 0;
		int no_of_rows_reached = 0;
		int TELE_Points_in_each_row_counter = 0;
		while(file >> row)
		{		
			if(row_index_only_in_TELE.size() ==0)
			{
				break;
			}

			if(i >= row_index_only_in_TELE.size())
			{	break;
			}

			if(no_of_rows_reached == row_index_only_in_TELE.at(i))
			{			
				if(transformed_TELE_coordinates_mapped_to_fusion_row.count(row_index_only_in_TELE[i]) == 0)
				{	
					i++;
					no_of_rows_reached++;
					continue;
				}
				itr = transformed_TELE_coordinates_mapped_to_fusion_row.find(row_index_only_in_TELE[i++]);
				std::ofstream myfile;
				myfile.open ("Output.csv" , std::ofstream::out | std::ofstream::app );
				if (!myfile.is_open())
				{
					std::cout << "Output operation not successfully performed\n";
					myfile.close();
					//return -1;
				}
				myfile <<	row[0] << ",";		
				myfile <<	row[1] << ",";		
				myfile <<	row[2] << ",";		
				myfile <<	row[3] << ",";		
				myfile <<	row[4] << ",";	
				myfile <<	(itr->second).left << ",";	
				myfile <<	(itr->second).bottom << ",";	
				myfile <<	(itr->second).right << ",";	
				myfile <<	(itr->second).top	<< ",";	
				myfile <<	(itr->second).xpos << ",";	
				myfile <<	(itr->second).ybottom	<< ",";	
				myfile <<	(itr->second).ytop << ",";	
				myfile <<	row[12]	<< ",";	
				myfile <<	row[13]	<< ",";	
				myfile <<	row[14]	<< ",";	
				myfile <<	row[15]	<< ",";	
				myfile <<	row[16]	<< ",";	
				myfile <<	row[17]	<< ",";	
				myfile <<	row[18]	<< ",";	
				myfile <<	row[19]	<< ",";	
				myfile <<	row[20] << ",";	
				myfile <<	row[21]	<< ",";	
				myfile <<	row[22]	<< ",";	
				myfile <<	row[23]	<< ",";	
				myfile <<	row[24]	<< ",";	
				myfile << std::endl;					
				myfile.close();						
			}
			no_of_rows_reached++;
		}

	}

}

void get_rows_from_csv(Camera type_of_Camera)
{
	if(type_of_Camera == FISHEYE)
	{
		std::ifstream file(FISHEYE_Csv);
		CSVRow row;
		int i = 0;
		int no_of_rows_reached = 0;
		int Fish_Eye_Points_in_each_row_counter = 0;
		while(file >> row)
		{	
			P.left = atof(row[5].c_str());
			P.bottom = atof(row[6].c_str());
			P.right	=	atof(row[7].c_str());
			P.top =	atof(row[8].c_str());
			P.xpos = atof(row[9].c_str());
			P.ybottom	=	atof(row[10].c_str());
			P.ytop = atof(row[11].c_str());
					
			fish_eye_coordinates_all_rows.insert(std::pair<double,Points>(i,P));
      i++;
		}
	}

	if(type_of_Camera == TELE)
	{
		std::ifstream file(TELE_Csv);
		CSVRow row;
		int i = 0;
		int no_of_rows_reached = 0;
		int TELE_Points_in_each_row_counter = 0;
		while(file >> row)
		{		
				P.left = atof(row[5].c_str());
				P.bottom	=	atof(row[6].c_str());
				P.right	=	atof(row[7].c_str());
				P.top =	atof(row[8].c_str());
				P.xpos =	atof(row[9].c_str());
				P.ybottom	=	atof(row[10].c_str());
				P.ytop = atof(row[11].c_str());

				TELE_coordinates_all_rows.insert(std::pair<double,Points>(i,P));
				i++;
		}

	}

}



void get_particular_csv_row_value(Camera type_of_Camera)
{
	if(type_of_Camera == FISHEYE)
	{
		std::ifstream file(FISHEYE_Csv);
		CSVRow row;
		int i = 0;
		int no_of_rows_reached = 0;
		int Fish_Eye_Points_in_each_row_counter = 0;
		while(file >> row)
		{	
			if(row_index_only_in_FISHEYE.size()==0)
			{
				break;
			}

			if(i >= row_index_only_in_FISHEYE.size())
			{	break;
			}			
			
			if(no_of_rows_reached == row_index_only_in_FISHEYE.at(i))
			{
				P.left = atof(row[5].c_str());
				P.bottom = atof(row[6].c_str());
				P.right =	atof(row[7].c_str());
				P.top	=	atof(row[8].c_str());
				P.xpos = atof(row[9].c_str());
				P.ybottom	=	atof(row[10].c_str());
				P.ytop = atof(row[11].c_str());

				Fish_Eye_Points_in_each_row_for_fusion.push_back(P);
					
				fish_eye_coordinates_mapped_to_fusion_row.insert(std::pair<double,Points>(
          row_index_only_in_FISHEYE[i],
          Fish_Eye_Points_in_each_row_for_fusion[Fish_Eye_Points_in_each_row_counter++])
          );

				i++;
			}
			no_of_rows_reached++;
		}
	}

	if(type_of_Camera == TELE)
	{
		std::ifstream file(TELE_Csv);
		CSVRow row;
		int i = 0;
		int no_of_rows_reached = 0;
		int TELE_Points_in_each_row_counter = 0;
		while(file >> row)
		{		
			if(row_index_only_in_TELE.size() ==0)
			{
				break;
			}

			if(i >= row_index_only_in_TELE.size())
			{	break;
			}

			if(no_of_rows_reached == row_index_only_in_TELE.at(i) )//throwing exception
			{
				P.left = atof(row[5].c_str());
				P.bottom = atof(row[6].c_str());
				P.right = atof(row[7].c_str());
				P.top =	atof(row[8].c_str());
				P.xpos = atof(row[9].c_str());
				P.ybottom =	atof(row[10].c_str());
				P.ytop = atof(row[11].c_str());

				TELE_Points_in_each_row_for_fusion.push_back(P);
				TELE_coordinates_mapped_to_fusion_row.insert(std::pair<double,Points>(row_index_only_in_TELE[i],
        TELE_Points_in_each_row_for_fusion[TELE_Points_in_each_row_counter++]));

				i++;
			}
			no_of_rows_reached++;
		}

	}

}







// Function for finding elements which are there
// in a[] but not in b[] . a - FISHEYE/TELE , b - main
//it contains the missing GI and row no. which is not in main but in other two Cam
void findMissing(std::vector<double>a,std::vector<double>b,Camera type_of_Camera)
{
		std::map<double, double>::iterator itr_main;
		std::map<double, double>::iterator itr_fe;
		std::map<double, double>::iterator itr_TELE;
		

	if(type_of_Camera == FISHEYE)
	{
		int counter_for_main_rows_remaining = 0;
		double value_to_insert=0;
		std::map<double,double>::iterator itr;
		std::map<double, double>::iterator itr1;
		std::map<double, double>::iterator itr_main_row;
		counter_for_main_rows_remaining = H_main.rowIndex;
		for (int i=0; i<a.size(); i++)
		{
			int j;
			for (j=0; j<b.size(); j++)
			{
				if (a[i] == b[j])
				{
					itr_main = main_xpos_mapped_to_grab_index.find(a[i]);
					itr_fe = fish_eye_xpos_mapped_to_grab_index.find(a[i]) ;
					if((itr_fe->second != 999 && itr_fe->second != 0) && (itr_main->second == 999 || itr_main->second == 0))
					{
						grab_index_only_in_a.push_back(itr_fe->first);
						itr = fish_eye_row_mapped_to_grab_index.find(itr_fe->first);
						row_index_only_in_FISHEYE.push_back(itr->second);
						
						value_to_insert = itr_fe->first;
	   				fish_eye_grab_index_mapped_to_row.insert(std::pair<double,double>(i,itr->first ));
		
						main_row_mapped_to_grab_index.erase(a[i]);
						counter_for_main_rows_remaining--;
					}
					break;
				}
			}
			if (j == b.size())
			{
						grab_index_only_in_a.push_back(a[i]);
						itr = fish_eye_row_mapped_to_grab_index.find(a[i]);
            //remember this row number is excel's row number- this is right- consistent
						row_index_only_in_FISHEYE.push_back(itr->second); 
						value_to_insert = a[i];
						fish_eye_grab_index_mapped_to_row.insert(std::pair<double,double>(itr->second,itr->first ));
			}
		}
		if(flag_for_CSV)
		{	
			get_particular_csv_row_value(FISHEYE);
		}
		
    }
	else if(type_of_Camera == TELE)
	{
		int counter_for_main_rows_remaining = 0;
		double value_to_insert=0;
		std::map<double, double>::iterator itr;
		for (int i=0; i<a.size(); i++)
		{
			int j;
			for (j=0; j<b.size(); j++)
			{
				if (a[i] == b[j])
				{
		
					itr_main = main_xpos_mapped_to_grab_index.find(a[i]);
					itr_TELE   = TELE_xpos_mapped_to_grab_index.find(a[i]);
					if((itr_TELE->second != 999 && itr_TELE->second != 0) && (itr_main->second == 999 || itr_main->second == 0))
					{
						grab_index_only_in_a.push_back(itr_TELE->first);
						itr = TELE_row_mapped_to_grab_index.find(itr_TELE->first);
						row_index_only_in_TELE.push_back(itr->second);
						value_to_insert = itr_TELE->first;
						TELE_grab_index_mapped_to_row.insert(std::pair<double,double>(i,itr->first));
						main_row_mapped_to_grab_index.erase(a[i]);
						counter_for_main_rows_remaining--;
					}
					break;
				}
			}
			if (j == b.size())
			{
			  grab_index_only_in_a.push_back(a[i]);
				itr = TELE_row_mapped_to_grab_index.find(a[i]);
				row_index_only_in_TELE.push_back(itr->second);
				
				value_to_insert = a[i];
				TELE_grab_index_mapped_to_row.insert(std::pair<double,double>(itr->second,itr->first));
				
				
			}
		}
		
		if(flag_for_CSV)
		{
				get_particular_csv_row_value(TELE);
		}
	
		
	}
	
}


int get_values_of_rows_in_sorted_main()
{
	int i=0;
	double value;
	int Counter1=0;
	int value_to_insert;
	image_coordinates_in_X_Y converted_points;

	std::ifstream       file;
	file.open (Main_Csv);
	if (!file.is_open())
	{
		std::cout << "Input file for Main Cam not present\n";
		file.close();
		return -1;
	}
	CSVRow row;
	Points P;
  while(file >> row)
    {		
		if(flag_for_VED == 0 && row.size() == NUMBER_OF_ATTRIBUTES_IN_VED)
		{	
			H_main.Keyframe	=	row[0];
			H_main.VCL = row[1];
			H_main.Perfect = row[2];
			H_main.camPort = row[3];
			H_main.originalID	=	row[4];
			H_main.left	=	row[5];
			H_main.bottom	=	row[6];
			H_main.right = row[7];
			H_main.top = row[8];
			H_main.xpos =	row[9];
			H_main.ybottom = row[10];
			H_main.ytop	=	row[11];
			H_main.type	=	row[12];
			H_main.truckType = row[13];
			H_main.lane	=	row[14];
			H_main.oncoming	=	row[15];
			H_main.CIPV	=	row[16];
			H_main.mobility	=	row[17];
			H_main.occluded	=	row[18];
			H_main.blinking	=	row[19];
			H_main.braking = row[20];
			H_main.grabIndex = row[21];
			H_main.enterHostLane = row[22];
			H_main.rearOcclusion = row[23];
			H_main.wheelVisibility = row[24];
			P.left = atoi(H_main.left.c_str());
			P.bottom = atoi(H_main.bottom.c_str());
			P.right = atoi(H_main.right.c_str());
			P.top = atoi(H_main.top.c_str());
			P.xpos = atoi(H_main.xpos.c_str());
			P.ybottom = atoi(H_main.ybottom.c_str());
			P.ytop = atoi(H_main.ytop.c_str());
		}
		else
		{
			std::cout<<"Wrong file for VED Main Cam is provided."<<std::endl;
			getchar();
			exit(0);
		}
		if(Counter1 == 0)
		{
			Global_Image_St_Row = atoi(H_main.Keyframe.c_str());
		}
		value_to_insert = atoi(H_main.grabIndex.c_str());
		sorted_main_coordinates_mapped_to_grab_index.insert(std::pair<double,Points>(value_to_insert,P));
		value = atoi(H_main.grabIndex.c_str());
		grab_index_sorted_main.push_back(value);
		Counter1++;
    }
	return 0;
}

double distance_Z_main ;
double Main_Min_Grab = 10000000;
double Main_Max_Grab = -1;


int get_values_of_rows_in_main()
{
	double value;
	int X1d, Y1d, X2d, Y2d, X3d, Y3d, X4d, Y4d;
	int X1u, Y1u, X2u, Y2u, X3u, Y3u, X4u, Y4u;
	char top_s[10], left_s[10], bottom_s[10], right_s[10], x_pos_s[10], y_top_s[10], y_bottom_s[10];

	double r_d = 0;
	double r_u = 0;
	int value_to_insert;
	std::string Inp_Left,Inp_Bot,Inp_Rig,Inp_Top,Inp_R1,Inp_R2,Inp_R3;
	image_coordinates_in_X_Y converted_points;

	double xpos_to_insert;

	std::ifstream       file;
	file.open (Main_Csv);
	if (!file.is_open())
	{
		std::cout << "Input file for Main Cam not present\n";
		file.close();
		return -1;
	}
	CSVRow              row;
	Points P;
	H_main.rowIndex = 0;
	while(file >> row)
    {		
		if(flag_for_VED == 0 && row.size() == NUMBER_OF_ATTRIBUTES_IN_VED)
		{	
			//Select a point in FISHEYE/TELE image (taken a bounding box bottom right)

			H_main.left	=	row[5];
			H_main.bottom	=	row[6];
			H_main.right = row[7];
			H_main.top = row[8];
			H_main.xpos	=	row[9];
			H_main.ybottom =	row[10];
			H_main.ytop =	row[11];

			 Inp_Left = row[5];
		   Inp_Bot = row[6];
			 Inp_Rig = row[7];
			 Inp_Top = row[8];
			 Inp_R1 = row[9];
			 Inp_R2 = row[10];
			 Inp_R3 = row[11];


			//Top Left Point
			X1d = atoi(row[5].c_str());
			Y1d = atoi(row[8].c_str());
			cartesian_to_image(STANDARD, converted_points, X1d, Y1d);
			X1d = converted_points.X;
			Y1d = converted_points.Y;
			r_d = calculate_r_d(STANDARD, X1d, Y1d);
			r_u = radial_distance_main(r_d);
			new_undistorted_points(STANDARD, r_u, r_d, X1d, Y1d, &X1u, &Y1u);
			// X1u += stCamParams.CamParams[STANDARD].OCX;
			// Y1u += stCamParams.CamParams[STANDARD].OCY;

			//Bottom Right Point
			X2d = atoi(row[7].c_str());
			Y2d = atoi(row[6].c_str());
			cartesian_to_image(STANDARD, converted_points, X2d, Y2d);
			X2d = converted_points.X;
			Y2d = converted_points.Y;
			r_d = calculate_r_d(STANDARD, X2d, Y2d);
			r_u = radial_distance_main(r_d);
			new_undistorted_points(STANDARD, r_u, r_d, X2d, Y2d, &X2u, &Y2u);
		
			// Xpos,ytop and Xpos,ybottom
			// Xpos,ytop Point
			X3d = atoi(row[9].c_str());
			Y3d = atoi(row[11].c_str());
			if(X3d != 999 && X3d != 0)
			{
				cartesian_to_image(STANDARD, converted_points, X3d, Y3d);
				X3d = converted_points.X;
				Y3d = converted_points.Y;
				r_d = calculate_r_d(STANDARD, X3d, Y3d);
				r_u = radial_distance_main(r_d);
				new_undistorted_points(STANDARD, r_u, r_d, X3d, Y3d, &X3u, &Y3u);
			}
				
	
			// Xpos,ybottom
			X4d = atoi(row[9].c_str());
			Y4d = atoi(row[10].c_str());
			if(X3d != 999 && X3d != 0)
			{
				cartesian_to_image(STANDARD, converted_points, X4d, Y4d);
				X4d = converted_points.X;
				Y4d = converted_points.Y;
				r_d = calculate_r_d(STANDARD, X4d, Y4d);
				r_u = radial_distance_main(r_d);
				new_undistorted_points(STANDARD, r_u, r_d, X4d, Y4d, &X4u, &Y4u);
			}
		

			image_coordinates_to_cartesian(STANDARD,converted_points,X1u,Y1u);
			X1u = converted_points.X;
			Y1u = converted_points.Y;
			image_coordinates_to_cartesian(STANDARD,converted_points,X2u,Y2u);
			X2u = converted_points.X;
			Y2u = converted_points.Y;

			if(X3d != 999 && X3d != 0)
			{
				image_coordinates_to_cartesian(STANDARD,converted_points,X3u,Y3u);
				X3u = converted_points.X;
				Y3u = converted_points.Y;
				image_coordinates_to_cartesian(STANDARD,converted_points,X4u,Y4u);
				X4u = converted_points.X;
				Y4u = converted_points.Y;
			}
	

			_itoa_s(Y1u, top_s, 10);
			_itoa_s(X1u, left_s, 10);
			_itoa_s(Y2u, bottom_s, 10);
			_itoa_s(X2u, right_s, 10);
			if(X3d != 999 && X3d != 0)
			{
				_itoa_s(X3u, x_pos_s   , 10);
				_itoa_s(Y3u, y_top_s   , 10);
				_itoa_s(Y4u, y_bottom_s, 10);
			}
      else{
				sprintf(x_pos_s, "999");
				sprintf(y_top_s, "999");
				sprintf(y_bottom_s, "999");
			}

			H_main.Keyframe	=	row[0];
			H_main.VCL = row[1];
			H_main.Perfect = row[2];
			H_main.camPort = row[3];
			H_main.originalID	=	row[4];
			H_main.left	=	std::string(left_s);
			H_main.bottom	=	std::string(bottom_s);
			H_main.right = std::string(right_s);
			H_main.top = std::string(top_s);
			H_main.xpos	=	std::string(x_pos_s);
			H_main.ybottom = std::string(y_bottom_s);
			H_main.ytop	=	std::string(y_top_s);
			H_main.type	=	row[12];
			H_main.truckType = row[13];
			H_main.lane	=	row[14];
			H_main.oncoming	=	row[15];
			H_main.CIPV	=	row[16];
			H_main.mobility	=	row[17];
			H_main.occluded	=	row[18];
			H_main.blinking	=	row[19];
			H_main.braking = row[20];
			H_main.grabIndex = row[21];
			H_main.enterHostLane = row[22];
			H_main.rearOcclusion = row[23];
			H_main.wheelVisibility = row[24];
			P.left = atoi(H_main.left.c_str());
			P.bottom = atoi(H_main.bottom.c_str());
			P.right = atoi(H_main.right.c_str());
			P.top = atoi(H_main.top.c_str());
			P.xpos = atoi(H_main.xpos.c_str());
			P.ybottom = atoi(H_main.ybottom.c_str());
			P.ytop = atoi(H_main.ytop.c_str());
		}
		else
		{
			std::cout<<"Wrong file for VED Main Cam is provided." <<std::endl;
			getchar();
			exit(0);
		}
		value_to_insert = atoi(H_main.grabIndex.c_str());
		if(value_to_insert > 0)
		{
			if(value_to_insert < Main_Min_Grab)
				Main_Min_Grab = value_to_insert;

		    if(value_to_insert > Main_Max_Grab)
				Main_Max_Grab = value_to_insert;

		}
		main_coordinates_mapped_to_grab_index.insert(std::pair<double,Points>(value_to_insert,P));
		main_row_mapped_to_grab_index.insert(std::pair<double,double>(value_to_insert,H_main.rowIndex));

		map_final_main_output_to_GrabIndex.insert(std::pair<double,VED_Header>(value_to_insert,H_main));
		map_final_main_output_to_row_number.insert(std::pair<double,VED_Header>(H_main.rowIndex,H_main));

		main_object_type_mapped_to_grab_index.insert(std::pair<int,int>(value_to_insert , atoi(H_main.type.c_str())));

		num_of_rows_main++;
		H_main.rowIndex++;
	
		H_main.left =	Inp_Left;
		H_main.bottom	=	Inp_Bot;
		H_main.right = Inp_Rig;
		H_main.top = Inp_Top;

		H_main.xpos =	Inp_R1;
		H_main.ybottom = Inp_R2;
		H_main.ytop	=	Inp_R3;

		map_main_output_to_GrabIndex.insert(std::pair<double,VED_Header>(value_to_insert,H_main));
		xpos_to_insert = atof(H_main.xpos.c_str());
		main_xpos_mapped_to_grab_index.insert(std::pair<double,double>(value_to_insert,xpos_to_insert));

		value = atoi(H_main.grabIndex.c_str());
		grab_index_main.push_back(value);
		
    }
	return 0 ;
}

double distance_Z_fish_eye;
double FE_Max_Grab = -1;
double FE_Min_Grab = 100000000;

int get_values_of_rows_in_fish_eye()
{
	double value;
	int value_to_insert;
	double xpos_to_insert;
	std::ifstream       file;
	file.open (FISHEYE_Csv);
	if (!file.is_open())
	{
		std::cout << "Input file for Fish Eye Cam not present\n";
		file.close();
		return -1;
	}
	CSVRow row;
	H_fish_eye.rowIndex = 0;
	while(file >> row)
    {	
		if(flag_for_VED == 0 && row.size() == NUMBER_OF_ATTRIBUTES_IN_VED)
		{	
			H_fish_eye.Keyframe	=	row[0];
			H_fish_eye.VCL = row[1];
			H_fish_eye.Perfect = row[2];
			H_fish_eye.camPort = row[3];
			H_fish_eye.originalID	=	row[4];
			H_fish_eye.left	=	row[5];
			H_fish_eye.bottom	=	row[6];
			H_fish_eye.right = row[7];
			H_fish_eye.top = row[8];
			H_fish_eye.xpos	=	row[9];
			H_fish_eye.ybottom = row[10];
			H_fish_eye.ytop	=	row[11];
			H_fish_eye.type	=	row[12];
			H_fish_eye.truckType	=	row[13];
			H_fish_eye.lane	=	row[14];
			H_fish_eye.oncoming	=	row[15];
			H_fish_eye.CIPV	=	row[16];
			H_fish_eye.mobility	=	row[17];
			H_fish_eye.occluded	=	row[18];
			H_fish_eye.blinking	=	row[19];
			H_fish_eye.braking	=	row[20];
			H_fish_eye.grabIndex	=	row[21];
			H_fish_eye.enterHostLane =	row[22];
			H_fish_eye.rearOcclusion = row[23];
			H_fish_eye.wheelVisibility =	row[24];
		}
		else
		{
			std::cout<<"Wrong file for VED Fish eye Cam is provided." <<std::endl;
			getchar();
			exit(0);
		}
	    value_to_insert = atoi(H_fish_eye.grabIndex.c_str());
		if(value_to_insert > 0)
		{
			if(value_to_insert < FE_Min_Grab)
				FE_Min_Grab = value_to_insert;

		    if(value_to_insert > FE_Max_Grab)
				FE_Max_Grab = value_to_insert;
		}

		fish_eye_row_mapped_to_grab_index.insert(std::pair<double,double>(value_to_insert,H_fish_eye.rowIndex));

		map_FISHEYE_output_to_GrabIndex.insert(
      std::pair<double,VED_Header>(value_to_insert,H_fish_eye)
      );

		xpos_to_insert = atof(H_fish_eye.xpos.c_str()) ;
		fish_eye_xpos_mapped_to_grab_index.insert(
      std::pair<double,double>(value_to_insert,xpos_to_insert)
      );

		
		fish_eye_object_type_mapped_to_grab_index.insert(
      std::pair<int,int>(value_to_insert,atoi(H_fish_eye.type.c_str()))
      );
	
		H_fish_eye.rowIndex++;
		num_of_rows_fish_eye++ ;
    value = atoi(H_fish_eye.grabIndex.c_str());
		grab_index_fish_eye.push_back(value) ;
		
    }
	return 0;
}

double distance_Z_TELE;
double TE_Min_Grab = 100000000;
double TE_Max_Grab = -1;

int get_values_of_rows_in_TELE()
{
	double value;
	int value_to_insert;
	double xpos_to_insert;
	std::ifstream file;
	file.open (TELE_Csv);
	if (!file.is_open())
	{
		std::cout << "Input file for TELE Cam not present\n";
		file.close();
		return -1;
	}
	CSVRow row;
	while(file >> row)
    {	
		if(flag_for_VED == 0 && row.size() == NUMBER_OF_ATTRIBUTES_IN_VED)
		{	
			H_TELE.Keyframe =	row[0];
			H_TELE.VCL = row[1];
			H_TELE.Perfect = row[2];
			H_TELE.camPort =	row[3];
			H_TELE.originalID	=	row[4];
			H_TELE.left	=	row[5];
			H_TELE.bottom	=	row[6];
			H_TELE.right = row[7];
			H_TELE.top = row[8];
			H_TELE.xpos	=	row[9];
			H_TELE.ybottom = row[10];
			H_TELE.ytop	=	row[11];
			H_TELE.type	=	row[12];
			H_TELE.truckType = row[13];
			H_TELE.lane	=	row[14];
			H_TELE.oncoming	=	row[15];
			H_TELE.CIPV	=	row[16];
			H_TELE.mobility	=	row[17];
			H_TELE.occluded	=	row[18];
			H_TELE.blinking	=	row[19];
			H_TELE.braking	=	row[20];
			H_TELE.grabIndex = row[21];
			H_TELE.enterHostLane =	row[22];
			H_TELE.rearOcclusion =	row[23];
			H_TELE.wheelVisibility =	row[24];
		}
		else
		{
			std::cout<<"Wrong file for VED TELE Cam is provided."<<std::endl;
			getchar();
			exit(0);
		}
		value_to_insert = atoi(H_TELE.grabIndex.c_str());
		if(value_to_insert > 0)
		{
			if(value_to_insert < TE_Min_Grab)
				TE_Min_Grab = value_to_insert;

		    if(value_to_insert > TE_Max_Grab)
				TE_Max_Grab = value_to_insert;
		}
  TELE_row_mapped_to_grab_index.insert(std::pair<double,double>(value_to_insert,H_TELE.rowIndex));

		map_TELE_output_to_GrabIndex.insert(std::pair<double,VED_Header>(value_to_insert,H_TELE));
			
		xpos_to_insert = atof(H_TELE.xpos.c_str());
		TELE_xpos_mapped_to_grab_index.insert(std::pair<double,double>(value_to_insert,xpos_to_insert));

		TELE_object_type_mapped_to_grab_index.insert(std::pair<int,int>(value_to_insert,atoi(H_TELE.type.c_str())));

		H_TELE.rowIndex++;
		num_of_rows_TELE++;
		value = atoi(H_TELE.grabIndex.c_str());
		grab_index_TELE.push_back(value);
		
    }
	return 0;
}

double DE_Min_Grab = 100000000;
double DE_Max_Grab = -1;
//Function for starting index of the image display
void Start_Row_Number_Set()
{

	if(TE_Min_Grab <= Main_Min_Grab)
	{
		Global_Start_Ind = Main_Min_Grab - TE_Min_Grab;
		return;
	}

	if((DE_Min_Grab <= Main_Min_Grab) && (TE_Min_Grab <= Main_Min_Grab))
	{
		Global_Start_Ind = DE_Min_Grab - TE_Min_Grab;
		return;
	}

}

double radial_distance_fish_eye(double r_d )
{
	double r_u;
	Cam = FISHEYE;
	double k1 = stCamParams.CamParams[Cam].distortionParams[0];
	double k2 = stCamParams.CamParams[Cam].distortionParams[1];
	double k3 = stCamParams.CamParams[Cam].distortionParams[2];
	double k4 = stCamParams.CamParams[Cam].distortionParams[3];
	double k5 = stCamParams.CamParams[Cam].distortionParams[4];
	double k6 = stCamParams.CamParams[Cam].distortionParams[5];
	double k7 = stCamParams.CamParams[Cam].distortionParams[6];
	double k8 = stCamParams.CamParams[Cam].distortionParams[7];
	double k9 = stCamParams.CamParams[Cam].distortionParams[8];

	double sK1 = (k1 >= 0) ? 1 : -1;
	double sK2 = (k2 >= 0) ? 1 : -1;
	double sK3 = (k3 >= 0) ? 1 : -1;
	double sK4 = (k4 >= 0) ? 1 : -1;
	double sK5 = (k5 >= 0) ? 1 : -1;
	double sK6 = (k6 >= 0) ? 1 : -1;
	double sK7 = (k7 >= 0) ? 1 : -1;
	double sK8 = (k8 >= 0) ? 1 : -1;
	double sK9 = (k9 >= 0) ? 1 : -1;

	double denominator = (
	1                  +
	sK1*pow(k1*r_d,2)  +
	sK2*pow(k2*r_d,4)  +
	sK3*pow(k3*r_d,6)  +
	sK4*pow(k4*r_d,8)  +
	sK5*pow(k5*r_d,10) +
	sK6*pow(k6*r_d,12) +
	sK7*pow(k7*r_d,14) +
	sK8*pow(k8*r_d,16) +
	sK9*pow(k9*r_d,18)
	);

	r_u = r_d/denominator;

	return r_u;
}

double radial_distance_main(double r_d )
{
	double r_u;
	Cam = STANDARD;
		
	double k1 = stCamParams.CamParams[Cam].distortionParams[0];
	double k2 = stCamParams.CamParams[Cam].distortionParams[1];
	double k3 = stCamParams.CamParams[Cam].distortionParams[2];

	double sK1 = (k1 >= 0) ? 1 : -1;
	double sK2 = (k2 >= 0) ? 1 : -1;
	double sK3 = (k3 >= 0) ? 1 : -1;
	
	double denominator = (
	1                  +
	sK1*pow(k1*r_d,2)  +
	sK2*pow(k2*r_d,4)  +
	sK3*pow(k3*r_d,6)  );

	r_u = r_d/denominator;

	return r_u;
}

double radial_distance_TELE(double r_d )
{
	double r_u;
	Cam = TELE;

	double k1 = stCamParams.CamParams[Cam].distortionParams[0];
	double k2 = stCamParams.CamParams[Cam].distortionParams[1];
	double sK1 = (k1 >= 0) ? 1 : -1;
	double sK2 = (k2 >= 0) ? 1 : -1;

	double denominator = (
	1                  +
	sK1*pow(k1*r_d,2)  +
	sK2*pow(k2*r_d,4) );

	r_u = r_d/denominator;

	return r_u;
}


// type_of_Camera- 0-STANDARD , 1- Fish-Eye , 2-TELE
double calculate_r_d(int type_of_Camera, int X, int Y)
{
	int new_X = X - stCamParams.CamParams[type_of_Camera].OCX;
	int new_Y = Y - stCamParams.CamParams[type_of_Camera].OCY;

	double r_d = sqrt((double)new_X*new_X + new_Y*new_Y);

	return r_d;
}

// type_of_Camera- 0-STANDARD , 1- Fish-Eye , 2-TELE
void new_undistorted_points(int type_of_Camera, double r_u, double r_d, double X_d, double Y_d,int *X_u ,int *Y_u)
{
	int image_height_Cam;
	if(type_of_Camera == TELE)
		 image_height_Cam = image_height_TELE;
	else
		 image_height_Cam = image_height;

	int new_X = X_d - stCamParams.CamParams[type_of_Camera].OCX;
	int new_Y = Y_d - stCamParams.CamParams[type_of_Camera].OCY;

	// The rotation and translation should be in cartesian coordinates at full resolution
	new_Y *= -1.0;

	// With respect to OCX and OCY
	int new_X_u = (r_u/r_d) * new_X;
	int new_Y_u = (r_u/r_d) * new_Y;

	new_Y_u *= -1.0;

	// Recalculate new_X_u, new_Y_u with respect to cx, cy
	*X_u = new_X_u + stCamParams.CamParams[type_of_Camera].OCX;
	*Y_u = new_Y_u + stCamParams.CamParams[type_of_Camera].OCY;
}

// Camera(TELE/F.E) specifying relative YPR  w.r.t. Main
void calculation_of_relative_YPR(Camera type_of_Camera, double *rel_alpha_roll ,double *rel_beta_yaw ,double *rel_gamma_pitch)
{
	double f_Cam = (type_of_Camera == FISHEYE) ? f_fe : f_TELE;

	int Cam_STANDARD = STANDARD;

	// TODO - Make it generic with respect to image height abd width of different Cameras
	int image_width_Cam  = image_width;
	int image_height_Cam = 0;
	if(type_of_Camera == FISHEYE )
		 image_height_Cam = image_height;
	else
		 image_height_Cam = image_height_TELE;

	int image_width_std = image_width;
	int image_height_std = image_height;

	int Cam_cy_1 = image_height_Cam - stCamParams.CamParams[type_of_Camera].cy;
	int difference_center_Cam_cy_1 = image_height_Cam/2 - Cam_cy_1;
	int std_cy_1 = image_height_std - stCamParams.CamParams[Cam_STANDARD].cy;
	int difference_center_std_cy_1 = image_height_std/2 - std_cy_1;

	int difference_center_Cam_cx = image_width_Cam/2 - stCamParams.CamParams[type_of_Camera].cx;
	int difference_center_std_cx = image_width_std/2 - stCamParams.CamParams[Cam_STANDARD].cx;
	
	// Calculation of rel_alpha_roll
	double alpha_roll = stCamParams.CamParams[type_of_Camera].RollAngle;
	double alpha_roll_std = stCamParams.CamParams[Cam_STANDARD].RollAngle;
	*rel_alpha_roll = alpha_roll - alpha_roll_std;

	double beta_yaw = atan((stCamParams.CamParams[type_of_Camera].Yaw + difference_center_Cam_cx) / f_Cam);
	double beta_yaw_std = atan((stCamParams.CamParams[Cam_STANDARD].Yaw + difference_center_std_cx) / f_std);
	*rel_beta_yaw = (beta_yaw - beta_yaw_std);

	double gamma_pitch = atan((stCamParams.CamParams[type_of_Camera].Horizon + difference_center_Cam_cy_1) / f_Cam);
	double gamma_pitch_std = atan((stCamParams.CamParams[Cam_STANDARD].Horizon + difference_center_std_cy_1) / f_std);
	*rel_gamma_pitch = (gamma_pitch - gamma_pitch_std);

}



void tranform_Cam_to_main(Camera type_of_Camera,
						  int x_Cam, int y_Cam,
						  int *final_Xm, int *final_Ym,double Z)
{

	double rel_gamma_pitch_Cam2m = 0.0f;
	double rel_beta_yaw_Cam2m = 0.0f;
	double rel_alpha_roll_Cam2m = 0.0f;

	// Cam to Main
	calculation_of_relative_YPR(type_of_Camera, &rel_alpha_roll_Cam2m, &rel_beta_yaw_Cam2m, &rel_gamma_pitch_Cam2m);

	cv::Mat Kfe = (cv::Mat_<double>(3, 3) <<
		f_fe, 0, stCamParams.CamParams[FISHEYE].cx,
		0, f_fe, stCamParams.CamParams[FISHEYE].cy,
		0, 0, 1);

	cv::Mat KTELE = (cv::Mat_<double>(3, 3) <<
		f_TELE, 0, stCamParams.CamParams[TELE].cx,
		0, f_TELE, stCamParams.CamParams[TELE].cy,
		0, 0, 1);

	cv::Mat Km = (cv::Mat_<double>(3, 3) <<
		f_std, 0, stCamParams.CamParams[STANDARD].cx,
		0, f_std, stCamParams.CamParams[STANDARD].cy,
		0, 0, 1);
	
	cv::Mat KCam = (type_of_Camera == FISHEYE) ? Kfe : KTELE;

	cv::Mat KmInv = Km.inv();
		
	cv::Mat RX1 = (cv::Mat_<double>(3, 3) <<
		1, 0, 0,
		0, cos(rel_gamma_pitch_Cam2m), -sin(rel_gamma_pitch_Cam2m),
		0, sin(rel_gamma_pitch_Cam2m), cos(rel_gamma_pitch_Cam2m));

	cv::Mat RY1 = (cv::Mat_<double>(3, 3) <<
		cos(rel_beta_yaw_Cam2m), 0, sin(rel_beta_yaw_Cam2m),
		0, 1, 0,
		-sin(rel_beta_yaw_Cam2m), 0, cos(rel_beta_yaw_Cam2m));

	cv::Mat RZ1 = (cv::Mat_<double>(3, 3) <<
		cos(rel_alpha_roll_Cam2m), -sin(rel_alpha_roll_Cam2m), 0,
		sin(rel_alpha_roll_Cam2m), cos(rel_alpha_roll_Cam2m), 0,
		0, 0, 1);

	// Camera axis are different from vehicle coordinate system
	cv::Mat R_all = RZ1 * RX1 * RY1;
	cv::Mat den = KCam * R_all * KmInv;

	/////////////////////////////////////////
	int X = 0;
	int Y = 0;
	double r_d = 0;
	double r_u = 0;
	int X_u = 0;
	int Y_u = 0;

	// Calculate in Camera Coordinate System
	double Cam_tz = stCamParams.CamParams[type_of_Camera].Tx;
	double Cam_tx = stCamParams.CamParams[type_of_Camera].Ty;
	double Cam_ty = stCamParams.CamParams[type_of_Camera].CamHt;

	double std_tz = stCamParams.CamParams[STANDARD].Tx;
	double std_tx = stCamParams.CamParams[STANDARD].Ty;
	double std_ty = stCamParams.CamParams[STANDARD].CamHt;

	cv::Mat tr = (cv::Mat_<double>(3, 1) <<
		Cam_tx - std_tx,
		Cam_ty - std_ty,
		Cam_tz - std_tz);

	// Select a point in FISHEYE/TELE image (taken a bounding box bottom right)
	X = x_Cam;
	Y = y_Cam;
	r_d = calculate_r_d(type_of_Camera, X, Y);

	if(type_of_Camera == FISHEYE)
	{
		r_u = radial_distance_fish_eye(r_d);
	}else
	{
		r_u = radial_distance_TELE(r_d);
	}
	
	X_u = 0;
	Y_u = 0;
	new_undistorted_points(type_of_Camera, r_u, r_d, X, Y, &X_u, &Y_u);
	cv::Mat mCam = (cv::Mat_<double>(3, 1) <<
		X_u,
		Y_u,
		1);
	float difference_in_ocy_and_center =  (stCamParams.CamParams[type_of_Camera].OCY - image_height/2);
	float Z_calculated= (stCamParams.CamParams[type_of_Camera].CamHt * stCamParams.CamParams[type_of_Camera].FocalLen)/ (Y_u + stCamParams.CamParams[type_of_Camera].Horizon + difference_in_ocy_and_center);
	Z = Z + Z_calculated ;
	
	cv::Mat num = mCam - ((KCam*tr) / Z);
	cv::Mat res = den.inv() * num;
	double val1 = res.at<double>(0);
	double val2 = res.at<double>(1);
	double val3 = res.at<double>(2);

	*final_Xm = val1/val3;
	*final_Ym = val2/val3;
}

int main(int argc, char **argv)
{


	int main_csv_file_present = -1;
	int fish_eye_csv_file_present	= -1;
	int TELE_csv_file_present	= -1;


	if(argc < 5)
	{
		std::cout<<"Please provide 0 for VED as first argument." <<std::endl;
		std::cout<<"Please provide 1 for CSV input files." <<std::endl;
		std::cout<<"Please provide name of main Cam csv file as third argument ,name of FISHEYE Cam csv file as fourth argument,";
		std::cout<<"name of TELE Cam csv file as fifth argument ,name of Camera calibration txt file as sixth argument."<<std::endl;
		std::cout<<"If any of the file not present please provide \"\""" in place of that argument"<<std::endl;
		std::cout<<"Press a key to exit."<<std::endl;
		getchar();
		return -1;
	}

  std::ofstream myfile;
	myfile.open("Output.csv", std::ofstream::out|std::ofstream::trunc);
	if(!myfile.is_open())
	{
		std::cout << "Output operation not successfully performed\n";
		myfile.close();
		return -1;
	}
	
	std::ofstream myfile1;
	myfile1.open("AllFused_Output.csv", std::ofstream::out|std::ofstream::trunc);
	if(!myfile1.is_open())
	{
		std::cout << "Output operation not successfully performed\n";
		myfile1.close();
		return -1;
	}

	if(argc>1)
	{
		 if (strcmp("0", argv[1]) == 0)
		 {
			 flag_for_VED = 0;
		 }
		 else
		 {
				std::cout<<"Please provide 0 for VED and 1 for PED as first argument"<<std::endl;
				std::cout<<"Press a key to exit."<<std::endl;
				getchar();
				return -1;
		 }
	}

	if(argc>2)
	{
		 if (strcmp("1", argv[2]) == 0)
		 {
			 flag_for_CSV = 1;
		 }
		 else
		 {
			  flag_for_CSV = 0;
		 }
	}


	if(flag_for_VED == 0)
	{
			if(flag_for_CSV)
			{				
				Main_Csv = argv[3];
				FISHEYE_Csv	= argv[4];	
				TELE_Csv = argv[5];
			}
			ConfigFile = argv[6];
	
			// Get Global Parameters
			get_cam_parameters(ConfigFile);
			f_fe = stCamParams.CamParams[FISHEYE].FocalLen;
			f_std = stCamParams.CamParams[STANDARD].FocalLen;
			f_TELE = stCamParams.CamParams[TELE].FocalLen;

			if(flag_for_CSV)
			{	main_csv_file_present = get_values_of_rows_in_main();
				fish_eye_csv_file_present = get_values_of_rows_in_fish_eye();
				TELE_csv_file_present = get_values_of_rows_in_TELE();
			}

			Start_Row_Number_Set(); // Global index for starting row
		
			if(main_csv_file_present == -1 && TELE_csv_file_present == -1 && fish_eye_csv_file_present == -1)
				return -1;
			if(fish_eye_csv_file_present != -1 || main_csv_file_present != -1)
				findMissing(grab_index_fish_eye , grab_index_main , FISHEYE);
			if(TELE_csv_file_present != -1 || main_csv_file_present != -1)
				findMissing(grab_index_TELE , grab_index_main , TELE);


      // write all the main output
			std::map<double,double>::iterator itr_search;
			std::map<double,VED_Header>::iterator itr_comapare;
			int i=0;
			for( i=0 ; i<map_final_main_output_to_row_number.size() ; i++)
			{
				if(map_final_main_output_to_row_number.count(i))
					itr_comapare = map_final_main_output_to_row_number.find(i);
				if(main_row_mapped_to_grab_index.count(grab_index_main[i]) > 0)
				{
					itr_search = main_row_mapped_to_grab_index.find(grab_index_main[i]);
					if(itr_search->second == itr_comapare->first)
						write_to_csv(itr_comapare->second);
				}
			}
			double Z_fe= 0;
			double Z_TELE= 0;
			point_in_image_coordinates f_e_final;
			point_in_image_coordinates TELE_final;
			std::map<double,double>::iterator itr;
			std::map<double,double>::iterator itr1;
			std::map<double,Points>::iterator itr_point;
			std::map<double,Points>::iterator itr_point1;
			
			int counter=0;
			int counter1=0;
			point_in_image_coordinates Point_f_e;
			Points P_fe;
			Points P_fe_transformed;

			point_in_image_coordinates Point_TELE;
			Points P_TELE;
			Points P_TELE_transformed;
			float absolute_translation_difference = 0;
			int grab_id = 0;
			if(fish_eye_csv_file_present != -1)
			while(counter < row_index_only_in_FISHEYE.size())
			{			
	
				itr = fish_eye_grab_index_mapped_to_row.find(row_index_only_in_FISHEYE[counter++]);
					grab_id  = itr->second;
					itr_point = fish_eye_coordinates_mapped_to_fusion_row.find(itr->first);

					P_fe.left = (itr_point->second).left;
					P_fe.bottom = (itr_point->second).bottom;
					P_fe.right = (itr_point->second).right;
					P_fe.top = (itr_point->second).top;
					P_fe.xpos = (itr_point->second).xpos;
					P_fe.ybottom = (itr_point->second).ybottom;
					P_fe.ytop = (itr_point->second).ytop;

					fish_eye_object_type_iterator = fish_eye_object_type_mapped_to_grab_index.find(grab_id);
					if(fish_eye_object_type_iterator->second == 1)
							absolute_translation_difference = 1.401;
					else
							absolute_translation_difference = 1.47;

					Z_fe = absolute_translation_difference;

					image_coordinates_in_X_Y converted_points;
					converted_points.X = 0;
					converted_points.Y = 0;
          // transforming cartesian coordinates (Top Left) to image coordinates
					cartesian_to_image(FISHEYE,converted_points, P_fe.left, P_fe.top);
					P_fe.left = converted_points.X;
					P_fe.top = converted_points.Y;
          // transforming cartesian coordinates (Bottom Right) to image coordinates
					cartesian_to_image(FISHEYE,converted_points, P_fe.right, P_fe.bottom);
					P_fe.right = converted_points.X;
					P_fe.bottom = converted_points.Y;
					if(P_fe.xpos != 999 && P_fe.xpos != 0)
					{
						cartesian_to_image(FISHEYE, converted_points, P_fe.xpos, P_fe.ytop);
						P_fe.xpos = converted_points.X;
						P_fe.ytop = converted_points.Y;
						cartesian_to_image(FISHEYE, converted_points, P_fe.xpos, P_fe.ybottom);
						P_fe.ybottom = converted_points.Y;
					}
					Point_f_e.bottomRight.X = P_fe.right;
					Point_f_e.bottomRight.Y = P_fe.bottom;
          Point_f_e.topLeft.X = P_fe.left;
					Point_f_e.topLeft.Y = P_fe.top;
          if(P_fe.xpos != 999 && P_fe.xpos != 0)
					{
						Point_f_e.xpos_ytop.X = P_fe.xpos;
						Point_f_e.xpos_ytop.Y = P_fe.ytop;

						Point_f_e.xpos_ybottom.X = P_fe.xpos;
						Point_f_e.xpos_ybottom.Y = P_fe.ybottom;
					}
					tranform_Cam_to_main(FISHEYE, Point_f_e.topLeft.X, Point_f_e.topLeft.Y, &f_e_final.topLeft.X, &f_e_final.topLeft.Y, Z_fe);
					tranform_Cam_to_main(FISHEYE, Point_f_e.bottomRight.X, Point_f_e.bottomRight.Y, &f_e_final.bottomRight.X, &f_e_final.bottomRight.Y, Z_fe);
					if(P_fe.xpos != 999 && P_fe.xpos != 0)
					{
						tranform_Cam_to_main(FISHEYE, Point_f_e.xpos_ytop.X, Point_f_e.xpos_ytop.Y, &f_e_final.xpos_ytop.X, &f_e_final.xpos_ytop.Y, Z_fe);
						tranform_Cam_to_main(FISHEYE, Point_f_e.xpos_ybottom.X, Point_f_e.xpos_ybottom.Y, &f_e_final.xpos_ybottom.X, &f_e_final.xpos_ybottom.Y, Z_fe);
					}
					image_coordinates_to_cartesian(STANDARD, converted_points, f_e_final.topLeft.X, f_e_final.topLeft.Y);
					P_fe_transformed.left = converted_points.X;
					P_fe_transformed.top = converted_points.Y;
					image_coordinates_to_cartesian(STANDARD, converted_points, f_e_final.bottomRight.X, f_e_final.bottomRight.Y);
					P_fe_transformed.right = converted_points.X;
					P_fe_transformed.bottom = converted_points.Y;
					if(P_fe.xpos != 999 && P_fe.xpos != 0)
					{
						image_coordinates_to_cartesian(STANDARD, converted_points, f_e_final.xpos_ytop.X, f_e_final.xpos_ytop.Y);
						P_fe_transformed.xpos = converted_points.X;
						P_fe_transformed.ytop = converted_points.Y;
						image_coordinates_to_cartesian(STANDARD, converted_points, f_e_final.xpos_ybottom.X, f_e_final.xpos_ybottom.Y);
						P_fe_transformed.xpos = converted_points.X;
						P_fe_transformed.ybottom = converted_points.Y;
					}
					if(P_fe.xpos == 999 || P_fe.xpos == 0)
					{
						P_fe_transformed.xpos = 999;
						P_fe_transformed.ytop = 999;
						P_fe_transformed.ybottom = 999;
					}
					transformed_fish_eye_coordinates_mapped_to_fusion_row.insert(std::pair<double,Points>(itr->first, P_fe_transformed));
			}
	
			if(TELE_csv_file_present!=-1)
			while(counter1 < row_index_only_in_TELE.size())
			{			
  			if(TELE_grab_index_mapped_to_row.count(row_index_only_in_TELE[counter1])> 0)
				{
					itr = TELE_grab_index_mapped_to_row.find(row_index_only_in_TELE[counter1++]);
					grab_id  = itr->second;
					itr_point = TELE_coordinates_mapped_to_fusion_row.find(itr->first);
          P_TELE.left =	(itr_point->second).left;
					P_TELE.bottom = (itr_point->second).bottom;
					P_TELE.right = (itr_point->second).right;
					P_TELE.top = (itr_point->second).top;
					P_TELE.xpos = (itr_point->second).xpos;
					P_TELE.ybottom = (itr_point->second).ybottom;
					P_TELE.ytop = (itr_point->second).ytop;
          TELE_object_type_iterator = TELE_object_type_mapped_to_grab_index.find(grab_id);
					if(TELE_object_type_iterator->second == 1)
							absolute_translation_difference = 1.401;
					else
							absolute_translation_difference = 1.47;
          Z_TELE = absolute_translation_difference;
          image_coordinates_in_X_Y converted_points;
					converted_points.X = 0;
					converted_points.Y = 0;
          // transforming cartesian coordinates (Top Left) to image coordinates
					cartesian_to_image(TELE,converted_points, P_TELE.left, P_TELE.top);
					P_TELE.left = converted_points.X;
					P_TELE.top = converted_points.Y;
          // transforming cartesian coordinates (Bottom Right) to image coordinates
					cartesian_to_image(TELE,converted_points, P_TELE.right, P_TELE.bottom);
					P_TELE.right = converted_points.X;
					P_TELE.bottom = converted_points.Y;
			    if(P_TELE.xpos != 999 && P_TELE.xpos != 0)
					{
						cartesian_to_image(TELE, converted_points, P_TELE.xpos, P_TELE.ytop);
						P_TELE.xpos = converted_points.X;
						P_TELE.ytop = converted_points.Y;
						cartesian_to_image(TELE, converted_points, P_TELE.xpos, P_TELE.ybottom);
						P_TELE.ybottom = converted_points.Y;
					}
					Point_TELE.bottomRight.X = P_TELE.right;
					Point_TELE.bottomRight.Y = P_TELE.bottom;
					Point_TELE.topLeft.X = P_TELE.left;
					Point_TELE.topLeft.Y = P_TELE.top;
					if(P_TELE.xpos != 999 && P_TELE.xpos != 0)
					{
						Point_TELE.xpos_ytop.X = P_TELE.xpos;
						Point_TELE.xpos_ytop.Y = P_TELE.ytop;
						Point_TELE.xpos_ybottom.X = P_TELE.xpos;
						Point_TELE.xpos_ybottom.Y = P_TELE.ybottom;
					}
          tranform_Cam_to_main(TELE, Point_TELE.topLeft.X, Point_TELE.topLeft.Y, &TELE_final.topLeft.X, &TELE_final.topLeft.Y, Z_TELE);
					tranform_Cam_to_main(TELE, Point_TELE.bottomRight.X, Point_TELE.bottomRight.Y, &TELE_final.bottomRight.X, &TELE_final.bottomRight.Y, Z_TELE);
					if(P_TELE.xpos != 999 && P_TELE.xpos != 0)
					{
						tranform_Cam_to_main(TELE, Point_TELE.xpos_ytop.X, Point_TELE.xpos_ytop.Y, &TELE_final.xpos_ytop.X, &TELE_final.xpos_ytop.Y, Z_TELE);
						tranform_Cam_to_main(TELE, Point_TELE.xpos_ybottom.X, Point_TELE.xpos_ybottom.Y, &TELE_final.xpos_ybottom.X, &TELE_final.xpos_ybottom.Y, Z_TELE);
					}
					image_coordinates_to_cartesian(STANDARD, converted_points, TELE_final.topLeft.X, TELE_final.topLeft.Y);
					P_TELE_transformed.left = converted_points.X;
					P_TELE_transformed.top = converted_points.Y;
					image_coordinates_to_cartesian(STANDARD, converted_points, TELE_final.bottomRight.X, TELE_final.bottomRight.Y);
					P_TELE_transformed.right = converted_points.X;
					P_TELE_transformed.bottom = converted_points.Y;
					if(P_TELE.xpos != 999 && P_TELE.xpos != 0)
					{
						image_coordinates_to_cartesian(STANDARD, converted_points, TELE_final.xpos_ytop.X, TELE_final.xpos_ytop.Y);
						P_TELE_transformed.xpos = converted_points.X;
						P_TELE_transformed.ytop = converted_points.Y;
						image_coordinates_to_cartesian(STANDARD, converted_points, TELE_final.xpos_ybottom.X, TELE_final.xpos_ybottom.Y);
						P_TELE_transformed.xpos = converted_points.X;
						P_TELE_transformed.ybottom = converted_points.Y;
					}
					if(P_TELE.xpos == 999 || P_TELE.xpos == 0)
					{
						P_TELE_transformed.xpos = 999;
						P_TELE_transformed.ytop = 999;
						P_TELE_transformed.ybottom = 999;
					}
					transformed_TELE_coordinates_mapped_to_fusion_row.insert(std::pair<double,Points>(itr->first, P_TELE_transformed));
				}
			}
	}
	
	if(fish_eye_csv_file_present != -1)
	{	
		write_to_particular_row(FISHEYE);
	}
	if(TELE_csv_file_present!=-1)
	{
		write_to_particular_row(TELE);
	}

    sorted_output();
	  Fused_Write();


//	getchar();
}

