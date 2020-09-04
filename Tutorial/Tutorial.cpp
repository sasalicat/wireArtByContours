// Tutorial.cpp : 定義主控台應用程式的進入點。
#include "stdafx.h"

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <ogldev_math_3d.h>
#include "OpenMesh\Core\IO\MeshIO.hh"
#include "OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh"
#include <stdlib.h>
#include <iostream>
#include <random>
//#include <glm>
//#include <time.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <stdio.h>
//#include <opencv2/core/core.hpp>
//#include <opencv2/opencv.hpp>
#include "script.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <fstream>
#include "polarMap_gdpc.h"
#include "tensorMap.h"
//#include "drawMap.h"
#include "imageConvert.h"
#include "intersectionOptimizer.h"
#include "geodesicPathFinder.h";
#include "IS_FeatureMatrix.h";
#include "meshCreater.h"
#include "WireComposition.h";
//#include <opencv2/core/core.hpp>
#include "GDPC\Mesh.h";
#include "GDPC\Generator.h"
//imgui

#include "imgui\imgui.h"
#include "imgui\imgui_impl_glut.h"
#include "imgui\imgui_impl_opengl2.h"

using namespace std;

//using namespace cv;
script *currentScript;

int SCR_WIDTH = 800;
int SCR_HEIGHT = 600;
//typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

struct Vertex
{
	Vector3f m_pos;
	Vector2f m_tex;
	Vector3f m_normal;
	
	Vertex() {}

	Vertex(Vector3f pos, Vector2f tex)
	{
		m_pos = pos;
		m_tex = tex;
		m_normal = Vector3f(0.0f, 0.0f, 0.0f);
	}
	Vertex(Vector3f pos, Vector2f tex,Vector3f norm)
	{
		m_pos = pos;
		m_tex = tex;
		m_normal = norm;
	}
};

static void AddShader(GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType)
{
	GLuint ShaderObj = glCreateShader(ShaderType);

	if (ShaderObj == 0) {
		fprintf(stderr, "Error creating shader type %d\n", ShaderType);
		//exit(0);
	}

	const GLchar* p[1];
	p[0] = pShaderText;
	GLint Lengths[1];
	Lengths[0] = strlen(pShaderText);
	glShaderSource(ShaderObj, 1, p, Lengths);
	glCompileShader(ShaderObj);
	GLint success;
	glGetShaderiv(ShaderObj, GL_COMPILE_STATUS, &success);
	if (!success) {
		GLchar InfoLog[1024];
		glGetShaderInfoLog(ShaderObj, 1024, NULL, InfoLog);
		fprintf(stderr, "Error compiling shader type %d: '%s'\n", ShaderType, InfoLog);
		//exit(1);
	}

	glAttachShader(ShaderProgram, ShaderObj);
}
static GLuint CompileShaders(char* pVSFileName, char* pFSFileName)
{
	GLuint ShaderProgram = glCreateProgram();

	if (ShaderProgram == 0) {
		fprintf(stderr, "Error creating shader program\n");
		//exit(1);
	}

	string vs, fs;

	if (!ReadFile(pVSFileName, vs)) {
		cout << "讀取" << pVSFileName << "失敗"<<endl;
		//exit(1);
	};

	if (!ReadFile(pFSFileName, fs)) {
		cout << "讀取" << pFSFileName << "失敗" << endl;
		//exit(1);
	};
	//cout <<"vs:"<< vs << endl;
	//cout << "fs:" << fs << endl;
	AddShader(ShaderProgram, vs.c_str(), GL_VERTEX_SHADER);
	AddShader(ShaderProgram, fs.c_str(), GL_FRAGMENT_SHADER);

	GLint Success = 0;
	GLchar ErrorLog[1024] = { 0 };

	glLinkProgram(ShaderProgram);
	glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &Success);
	if (Success == 0) {
		glGetProgramInfoLog(ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
		fprintf(stderr, "Error linking shader program: '%s'\n", ErrorLog);
		//exit(1);
	}

	glValidateProgram(ShaderProgram);
	glGetProgramiv(ShaderProgram, GL_VALIDATE_STATUS, &Success);
	if (!Success) {
		glGetProgramInfoLog(ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
		fprintf(stderr, "Invalid shader program: '%s'\n", ErrorLog);
		//exit(1);
	}
	cout << "編譯完成"<<endl;
	return ShaderProgram;
}

class showModel:public script
{	
public:
	GLuint VBO[2];
	void onInit(){
		cout << "showModel onInit"<<endl;
		//Vector3f Vertices[3];
		//Vertices[0] = Vector3f(-1.0f, -1.0f, 0.0f);
		//Vertices[1] = Vector3f(1.0f, -1.0f, 0.0f);
		//Vertices[2] = Vector3f(0.0f, 1.0f, 0.0f);
		//mesh.request_vertex_normals();
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		const char* filename = "../Resources/bunny.obj";
		if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
		{
			cerr << "Error: Cannot read mesh from " << filename << endl;
			return;
		}
		else
		{
			cout << "load mesh成功" << endl;
		}
		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}
		mesh.update_normals();
		/*
		if (!mesh.has_vertex_normals()){
			std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
			mesh.normal
			//return;
		}
		else{
			cout << "normal 存在" << endl;
		}*/
		//mesh.calc_vertex_normal();

		mesh.request_vertex_texcoords2D();
		//vector<GLfloat> Vertices = vector<GLfloat>(mesh.n_faces()*9);
		
		Vector3f *Vertices = new Vector3f[mesh.n_faces() * 3];
		Vector3f *Normals = new Vector3f[mesh.n_faces() * 3];
		for (int f = 0; f < mesh.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				Vertices[f*3+vindex] = Vector3f(pos[0], pos[1], pos[2]);
				MyMesh::Normal norm = mesh.normal(*fv_it);
				Normals[f * 3 + vindex] = Vector3f(norm[0],norm[1],norm[2]);
				//cout << "面" << f << "-norm" << vindex << ":(" << norm[0] << " ," << norm[1] << " ," << norm[2] << ")" << endl;

				//cout << "面" << f << "-point" << vindex << ":(" << pos[0] << " ," << pos[1] << " ," << pos[2] << ")" << endl;

				//float* pos = mesh.point(*fv_it).data();


				//Vertices[f * 9 + vindex * 3] = pos[0];
				//Vertices[f * 9 + vindex * 3 + 1] = pos[1];
				//Vertices[f * 9 + vindex * 3 + 2] = pos[2];
				vindex++;
			}
			//cout << endl;
		}
		/*
		Vector3f* Vertices = new Vector3f[6];
		Vertices[0] = Vector3f(-1.0f, -1.0f, .0f);
		Vertices[1] = Vector3f(.0f, -1.0f, .0f);
		Vertices[2] = Vector3f(.0f, 1.0f, .0f);
		Vertices[3] = Vector3f(1.0f, -1.0f, .0f);
		Vertices[4] = Vector3f(.0f, -1.0f, .0f);
		Vertices[5] = Vector3f(.0f, 1.0f, .0f);*/
		//cout << "Vertices size 为:" << sizeof(Vertices);

		/*glGenTextures(1, &frameTexHandle);
		glBindTexture(GL_TEXTURE_2D, frameTexHandle);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, frameTexHandle, 0);
		*/



		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Normals[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

		//深度buffer
		/*
		glGenFramebuffers(1, &fbHandle);
		glBindFramebuffer(GL_FRAMEBUFFER, fbHandle);
		glGenTextures(1, &textureID);
		glActiveTexture(GL_TEXTURE0 + 100);
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, SCR_WIDTH, SCR_WIDTH, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, textureID, 0);
		*/
		
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		//glEnable(GL_DEPTH_TEST);
		shaderHandle = CompileShaders("lighting.vs", "dir_light.fs");
		glUseProgram(shaderHandle);
		
		Vector3f vertex[6];
		vertex[0] = Vector3f(-1.0, -1.0, .0);
		vertex[1] = Vector3f(-1.0, 1.0, .0);
		vertex[2] = Vector3f(1.0, -1.0, .0);
		vertex[3] = Vector3f(-1.0, 1.0, .0);
		vertex[4] = Vector3f(1.0, -1.0, .0);
		vertex[5] = Vector3f(1.0, 1.0, .0);

		Vector2f coord[6];
		coord[0] = Vector2f(.0, 1.0);
		coord[1] = Vector2f(.0, .0);
		coord[2] = Vector2f(1.0, 1.0);
		coord[3] = Vector2f(.0, .0);
		coord[4] = Vector2f(1.0, 1.0);
		coord[5] = Vector2f(1.0, .0);
		//第二个pass的资料设置
		GLuint VBO_pass2[2];
		glGenVertexArrays(1, &VAO_pass2);
		glBindVertexArray(VAO_pass2);

		glGenBuffers(2, VBO_pass2);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, 6 * sizeof(Vector3f), &vertex[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[1]);
		glBufferData(GL_ARRAY_BUFFER, 6 * sizeof(Vector2f), &coord[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[1]);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
		
		//创建图片
		//深度buffer
		/*glGenTextures(1, &frameDepthHandle);
		glActiveTexture(GL_TEXTURE0 + 99);
		glBindTexture(GL_TEXTURE_2D, frameDepthHandle);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, SCR_WIDTH, SCR_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, frameDepthHandle, 0);
		*/

		//颜色buffer
		
		glGenTextures(1, &textureID);

		// "Bind" the newly created texture : all future texture functions will modify this texture
		//glActiveTexture(GL_TEXTURE0+100);
		glBindTexture(GL_TEXTURE_2D, textureID);
		//glActiveTexture(GL_TEXTURE0 + 100);
		// Give the image to OpenGL
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		//创建framebuffer
		
		glGenFramebuffers(1, &fbHandle);		
		glBindFramebuffer(GL_FRAMEBUFFER, fbHandle);		
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);

		
		//GLuint depthTextureID;
		glGenTextures(1, &depthID);
		glBindTexture(GL_TEXTURE_2D, depthID);
		
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		//glTexStorage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, SCR_WIDTH, SCR_HEIGHT);
		//glTexStorage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32F, SCR_WIDTH, SCR_HEIGHT);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, SCR_WIDTH, SCR_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthID, 0);

		//GLenum bufs[] = {GL_COLOR_ATTACHMENT0};
	//	glDrawBuffers(1, bufs);

		cerr << glCheckFramebufferStatus(GL_FRAMEBUFFER);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		{
			cout << "frame buffer grejkgojeojop!" << endl;
		}


		glBindTexture(GL_TEXTURE_2D, textureID);
		//unsigned int rbo;
		//glGenRenderbuffers(1, &rbo);
		//glBindRenderbuffer(GL_RENDERBUFFER, rbo);
		//glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, SCR_WIDTH, SCR_HEIGHT);
		//glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);
		//glBindRenderbuffer(GL_RENDERBUFFER, 0);
		
		
		pass2_shader = CompileShaders("tex.vs", "tex.fs");

		//glUseProgram(pass2_shader);
		cout << "准备完成" << endl;
	}
	void onRender(){
		cout << "render...";
		float speed = 0.01;

		//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		glBindFramebuffer(GL_FRAMEBUFFER, fbHandle);
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glBindVertexArray(VAO);
		glUseProgram(shaderHandle);
		//glUseProgram(shaderHandle);
		//angle += speed;
		offset -= speed;
		glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0, 1, 0));
		glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, offset));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 0, 0, z }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		GLuint loc= glGetUniformLocation(shaderHandle,"gWVP");
		assert(loc != 0xFFFFFFFF);
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(viewMat*translate));
		
		glm::vec3 lightDir = glm::vec3(1, 0, 0);
		GLuint loc_light = glGetUniformLocation(shaderHandle, "light_direction");
		assert(loc_light != 0xFFFFFFFF);
		glUniform3f(loc_light, lightDir[0], lightDir[1], lightDir[2]);
		//glUniform4f()
		//glClear(GL_COLOR_BUFFER_BIT);

		glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces()*9);
		glBindVertexArray(0);
		//glDisable(GL_DEPTH_TEST);
		//glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);


		// Clear all relevant buffers
		//glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set clear color to white (not really necessery actually, since we won't be able to see behind the quad anyways)
		//glClear(GL_COLOR_BUFFER_BIT);
		//glDisable(GL_DEPTH_TEST);
		glBindFramebuffer(GL_FRAMEBUFFER,0);
		glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
		glDisable(GL_DEPTH_TEST);
		glClear(GL_COLOR_BUFFER_BIT);
		glUseProgram(pass2_shader);
		glBindVertexArray(VAO_pass2);


		//glBindFramebuffer(GL_FRAMEBUFFER,fbHandle);
		//glBindTexture(GL_TEXTURE_2D, textureID);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, depthID);
		GLuint loc_tex = glGetUniformLocation(pass2_shader, "ourTexture");
		//assert(loc_tex != 0xFFFFFFFF);
		glUniform1i(loc_tex, 0);
		GLuint loc_near = glGetUniformLocation(pass2_shader,"zNear");
		glUniform1f(loc_near, nearZ);
		GLuint loc_far = glGetUniformLocation(pass2_shader, "zFar");
		glUniform1f(loc_far, farZ);

		glDrawArrays(GL_TRIANGLES,0,6);
		//glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		glBindVertexArray(0);
		//z += speed;
		//glEnableVertexAttribArray(0);
		//glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

private:
	MyMesh mesh;
	GLuint shaderHandle;
	GLuint pass2_shader;
	GLuint frameTexHandle;
	GLuint frameDepthHandle;
	float z = 1;
	float offset = 0;
	float angle = 0;
	float nearZ = 0.1f;
	float farZ = 30.0f;
	float threshold = 0.35f;
	GLuint VAO;
	GLuint VAO_pass2;
	GLuint fbHandle;
	GLuint textureID;
	GLuint depthID;
};

class showTexture:public script{
public:
	void onInit(){
		Vector3f vertex[6];
		vertex[0] = Vector3f(-1.0,-1.0, .0);
		vertex[1] = Vector3f(-1.0, 1.0, .0);
		vertex[2] = Vector3f(1.0, -1.0, .0);
		vertex[3] = Vector3f(-1.0, 1.0, .0);
		vertex[4] = Vector3f(1.0, -1.0, .0);
		vertex[5] = Vector3f(1.0, 1.0, .0);

		Vector2f coord[6];
		coord[0] = Vector2f(.0, 1.0);
		coord[1] = Vector2f(.0, .0);
		coord[2] = Vector2f(1.0, 1.0);
		coord[3] = Vector2f(.0, .0);
		coord[4] = Vector2f(1.0, 1.0);
		coord[5] = Vector2f(1.0, .0);
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, 6* sizeof(Vector3f), &vertex[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, 6* sizeof(Vector2f), &coord[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
		int width, height, nrChannels;
		unsigned char *data = stbi_load("../Resources/texture.jpg", &width, &height, &nrChannels, 0);
		GLuint textureID;
		glGenTextures(1, &textureID);

		// "Bind" the newly created texture : all future texture functions will modify this texture
		glBindTexture(GL_TEXTURE_2D, textureID);

		// Give the image to OpenGL
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		shaderHandle = CompileShaders("tex.vs", "tex.fs");
		glUseProgram(shaderHandle);
		cout << "onInit 完成"<<endl;
	}
	void  onRender(){
		glBindVertexArray(VAO);
		glClear(GL_COLOR_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);
	}
private:
	GLuint shaderHandle;
	GLuint VAO;
	GLuint VBO[2];
};

class  showSM:public script
{
public:
	//GLfloat* vposArray;
	//GLfloat* vcoordArray;
	//GLfloat* vnormArray;
	GLuint vaoHandle;
	GLint shaderHandle;
	GLuint vboHandles[3];
	void onInit(){

		shaderHandle = CompileShaders("lighting.vs", "showNormal.fs");
		glUseProgram(shaderHandle);
		OpenMesh::IO::Options opt;
		const char* filename = "../Resources/bunny.obj";
		if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
		{
			cerr << "Error: Cannot read mesh from " << filename << endl;
			return;
		}
		else
		{
			cout << "load mesh成功"<<endl;
		}
		mesh.request_vertex_normals();
		if (!mesh.has_vertex_normals()){
			std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
			return;
		}
		else{
			cout << "normal 存在" << endl;
		}
		mesh.request_vertex_texcoords2D();
		/*if (!mesh.has_halfedge_texcoords2D())
		{
			std::cerr << "no texcoord2D exist!\n";
			return;
		}
		else
		{
			cout << "貼圖坐標存在" << endl;
		}
		*/
		cout << "模型面數:" << mesh.n_faces() << endl;
		//int face_num = mesh.n_faces();
		//int pos_num = face_num * 9;
		//從模型里加載頂點信息
		vector<GLfloat>  locals = vector<GLfloat>(mesh.n_faces() * 9);
		vector<GLfloat> coords = vector<GLfloat>(mesh.n_faces() * 6);
		vector<GLfloat> norms = vector<GLfloat>(mesh.n_faces() * 9);
		cout << "locals  init size:" << locals.size()<<endl;
		for (int f = 0; f < mesh.n_faces(); f++){
			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh);fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				

				locals[f*9 + vindex*3] =pos[0];
				locals[f * 9 + vindex * 3+1] = pos[1];
				locals[f * 9 + vindex * 3+2] = pos[2];
				//cout << "locals:" << pos[0] << "," << pos[1]<<","<<pos[2]<<"  ";
				//OpenMesh::Vec2f coord = mesh.texcoord2D(*fv_it);
				//coords.push_back(coord[0]);
				//coords.push_back(coord[1]);
				coords[f * 6 + vindex * 2] = 0;
				coords[f * 6 + vindex * 2 + 1] = 0;

				MyMesh::Normal norm= mesh.normal(*fv_it);
				norms[f * 9 + vindex * 3]=(norm.data()[0]);
				norms[f * 9 + vindex * 3+1]=(norm.data()[1]);
				norms[f * 9 + vindex * 3+2]=(norm.data()[2]);
				//Vertex vpos = Vertex(Vector3f(pos[0],pos[1],pos[2]),Vector2f(coord[0],coord[1]));
				vindex++;
			}
		}//加載完成

		//GLuint vboHandles[3];
		glGenBuffers(3, vboHandles);
		cout << "test begin" << endl;
		cout << "local size:" + locals.size() << " 前30個數字為:"<<endl;
		GLfloat* localArray = &locals[0];
		for (int i = 0; i < 30; i++){
			cout << localArray[i] << " ,";
		}
		cout << "local size:" << locals.size()*sizeof(GLfloat) << endl;
		posBufferHandle = vboHandles[0];
		glBindBuffer(GL_ARRAY_BUFFER, posBufferHandle);
		glBufferData(GL_ARRAY_BUFFER, locals.size()*sizeof(GLfloat), &locals[0], GL_STATIC_DRAW);
		
		texBufferHandle = vboHandles[1];
		glBindBuffer(GL_ARRAY_BUFFER, texBufferHandle);
		glBufferData(GL_ARRAY_BUFFER, coords.size()*sizeof(GLfloat), &coords[0], GL_STATIC_DRAW);

		normBufferHandle = vboHandles[2];
		glBindBuffer(GL_ARRAY_BUFFER, normBufferHandle);
		glBufferData(GL_ARRAY_BUFFER, norms.size()*sizeof(GLfloat), &norms[0], GL_STATIC_DRAW);
		glUseProgram(shaderHandle);
		glGenVertexArrays(1, &vaoHandle);
		glBindVertexArray(vaoHandle);

		cout << "绑定VAO:" << vaoHandle;
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);

		glBindBuffer(GL_ARRAY_BUFFER, posBufferHandle);

		glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(GLubyte*)NULL);

		glBindBuffer(GL_ARRAY_BUFFER, texBufferHandle);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (GLubyte*)NULL);

		glBindBuffer(GL_ARRAY_BUFFER, normBufferHandle);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte*)NULL); 
		glUseProgram(shaderHandle);
		
		

	}

	void  onRender(){
		cout << "使用VAO:" << vaoHandle<<endl;
		//cout << "onrender..." << endl;
		glClear(GL_COLOR_BUFFER_BIT);


		glBindVertexArray(vaoHandle);

		glUseProgram(shaderHandle);
		glBindBuffer(GL_ARRAY_BUFFER, vboHandles[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte*)NULL);
		glEnableVertexAttribArray(0);
		/*glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, posBufferHandle);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte*)NULL);

		glBindBuffer(GL_ARRAY_BUFFER, texBufferHandle);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (GLubyte*)NULL);

		glBindBuffer(GL_ARRAY_BUFFER, normBufferHandle);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte*)NULL);*/
		//glUseProgram(shaderHandle);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);

		//glutSwapBuffers();
	}
private: 
	MyMesh mesh;
	GLuint posBufferHandle;
	GLuint texBufferHandle;
	GLuint normBufferHandle;
};


class tt2 :public script{
public:
	GLuint VBO;
	void onInit(){
		string sphere_path = "../Resources/uvSphere_16.obj";
		Vector3f Vertices[1];
		Vertices[0] = Vector3f(.0f, .0f, .0f);
		//Vertices[1] = Vector3f()
		int center = 296;
		typedef DGPC::Vector3<double> Point;
		typedef DGPC::MeshOM<Point> Mesh;
		typedef DGPC::Generator<Mesh> DGPCgenerator;
		//typedef DGPC::Generator<MyMesh> DGPCgenerator_2;
		Mesh mymesh;
		mymesh.openOBJ(sphere_path.c_str());
		//MyMesh mesh2;
		//DGPCgenerator_2 dgpc(mesh2);
		DGPCgenerator dgpc(mymesh);
		dgpc.setNodeSource(center);
		int last_node = dgpc.run();
		cout << "DGPC 完成!!!!!!" << endl;
		std::cout << "Computed distances until node " << last_node << std::endl;
		std::cout << std::endl;

		std::cout << "i      r      theta" << std::endl;
		std::cout << "-------------------" << std::endl;
		for (int i = 0; i < mymesh.n_vertices(); i++) {
			const double r = dgpc.getDistance(i);
			if (r < std::numeric_limits<double>::max()) {
				const double theta = dgpc.getAngle(i);
				std::cout << i << "    " << r << "    " << theta << std::endl;
			}
		}
		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
	}
	void onRender(){
		cout << "tt2 onrender";
		glClear(GL_COLOR_BUFFER_BIT);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glDrawArrays(GL_POINTS, 0, 1);

		glEnableVertexAttribArray(0);
	}
};
class tt3 :public script{
public:
	GLuint VBO;
	void onInit(){
		Vector3f Vertices[3];
		Vertices[0] = Vector3f(-1.0f, -1.0f, .0f);
		Vertices[1] = Vector3f(1.0f, -1.0f, .0f);
		Vertices[2] = Vector3f(.0f,1.0f,.0f);
		//Vertices[1] = Vector3f()

		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
	}
	void onRender(){
		glClear(GL_COLOR_BUFFER_BIT);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		//glDrawBuffer(GL_DEPTH);
		//glReadBuffer(GL_COLOR);
		glDrawArrays(GL_TRIANGLES, 0, 3);

		glEnableVertexAttribArray(0);
	}
};

class tt4 :public script{
	GLuint VBO;
	void onInit(){
		
		Vector3f Vertices[3];
		Vertices[0] = Vector3f(-1.0f, -1.0f, 0.0f);
		Vertices[1] = Vector3f(1.0f, -1.0f, 0.0f);
		Vertices[2] = Vector3f(0.0f, 1.0f, 0.0f);

		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
		glUseProgram(CompileShaders("shader.vs", "shader.fs"));
	}
	void onRender(){
		glClear(GL_COLOR_BUFFER_BIT);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glDrawArrays(GL_TRIANGLES, 0, 3);

		glEnableVertexAttribArray(0);
	}

};
class tt5 :public script{
	GLuint VBO;
	GLuint gScaleLocation;
	float scale = 0.f;
	void onInit(){

		Vector3f Vertices[3];
		Vertices[0] = Vector3f(-1.0f, -1.0f, 0.0f);
		Vertices[1] = Vector3f(1.0f, -1.0f, 0.0f);
		Vertices[2] = Vector3f(0.0f, 1.0f, 0.0f);

		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
		GLuint ShaderProgram = CompileShaders("shader2.vs", "shader2.fs");
		glUseProgram(ShaderProgram);
		gScaleLocation = glGetUniformLocation(ShaderProgram, "gScale");
		assert(gScaleLocation != 0xFFFFFFFF);
	}
	void onRender(){
		glClear(GL_COLOR_BUFFER_BIT);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glDrawArrays(GL_TRIANGLES, 0, 3);

		glEnableVertexAttribArray(0);
		scale += 0.01f;
		glUniform1f(gScaleLocation, scale);
	}
};
class  deferredShading:public script
{
	void onInit(){
		GLuint gBuffer;
		glGenFramebuffers(1, &gBuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);
		GLuint gPosition, gNormal, gColorSpec;

		// - 位置颜色缓冲
		glGenTextures(1, &gPosition);
		glBindTexture(GL_TEXTURE_2D, gPosition);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gPosition, 0);

			// - 法线颜色缓冲
		glGenTextures(1, &gNormal);
		glBindTexture(GL_TEXTURE_2D, gNormal);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gNormal, 0);

		// - 颜色 + 镜面颜色缓冲
		glGenTextures(1, &gColorSpec);
		glBindTexture(GL_TEXTURE_2D, gColorSpec);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, gColorSpec, 0);

		// - 告诉OpenGL我们将要使用(帧缓冲的)哪种颜色附件来进行渲染
		GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
		glDrawBuffers(3, attachments);
	}
};
class tt7 :public script{
	GLuint VBO;
	GLuint gScaleLocation;
	float scale = 0.f;
	void onInit(){

		Vector3f Vertices[3];
		Vertices[0] = Vector3f(-1.0f, -1.0f, 0.0f);
		Vertices[1] = Vector3f(1.0f, -1.0f, 0.0f);
		Vertices[2] = Vector3f(0.0f, 1.0f, 0.0f);

		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
		GLuint ShaderProgram = CompileShaders("shader3.vs", "shader3.fs");
		glUseProgram(ShaderProgram);
		gScaleLocation = glGetUniformLocation(ShaderProgram, "gWorld");
		assert(gScaleLocation != 0xFFFFFFFF);
	}
	void onRender(){
		glClear(GL_COLOR_BUFFER_BIT);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glDrawArrays(GL_TRIANGLES, 0, 3);

		glEnableVertexAttribArray(0);
		float translate[4][4];
		translate[0][0] = 1.0f;
		translate[1][0] = 0.0f;
		translate[2][0] = 0.0f;
		translate[3][0] = 0.0f;
		translate[0][1] = .0f;
		translate[1][1] = 1.0f;
		translate[2][1] = 0.0f;
		translate[3][1] = 0.0f;
		translate[0][2] = .0f;
		translate[1][2] = 0.0f;
		translate[2][2] = 1.0f;
		translate[3][2] = 0.0f;
		translate[0][3] = 0.0f;
		translate[1][3] = 0.0f;
		translate[2][3] = 1.0f;
		translate[3][3] = 0.0f;
		glUniform1f(gScaleLocation, scale);
	}

};
MyMesh normalizeMesh(MyMesh input){
	MyMesh::Point meshCenter = MyMesh::Point(.0, .0, .0);

	//遍历所有模型顶点加总获得平均值
	for (int i = 0; i < input.n_vertices(); i++){
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		MyMesh::Point point = input.point(vh);
		meshCenter += point;
	}
	meshCenter /= input.n_vertices();
	cout << "avg  mesh center:(" << meshCenter[0] << "," << meshCenter[1] << "," << meshCenter[2] << ")"<<endl;
	//更新所有的顶点到相对于中心点的位置
	float maxLength = 0;
	for (int i = 0; i < input.n_vertices(); i++){
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		MyMesh::Point point = input.point(vh);
		float len = (point - meshCenter).length();
		if (len > maxLength)
		{
			maxLength = len;
		}
	}
	for (int i = 0; i < input.n_vertices(); i++){
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		MyMesh::Point point = input.point(vh);
		input.set_point(vh, (point - meshCenter) /maxLength);
	}
	
	input.update_normals();
	return input;
}
MyMesh normalizeMesh(MyMesh input,float scale){
	MyMesh::Point meshCenter = MyMesh::Point(.0, .0, .0);

	//遍历所有模型顶点加总获得平均值
	for (int i = 0; i < input.n_vertices(); i++){
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		MyMesh::Point point = input.point(vh);
		meshCenter += point;
	}
	meshCenter /= input.n_vertices();
	cout << "avg  mesh center:(" << meshCenter[0] << "," << meshCenter[1] << "," << meshCenter[2] << ")" << endl;
	//更新所有的顶点到相对于中心点的位置
	float maxLength = 0;
	for (int i = 0; i < input.n_vertices(); i++){
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		MyMesh::Point point = input.point(vh);
		float len = (point - meshCenter).length();
		if (len > maxLength)
		{
			maxLength = len;
		}
	}
	for (int i = 0; i < input.n_vertices(); i++){
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		MyMesh::Point point = input.point(vh);
		input.set_point(vh, ((point - meshCenter) / maxLength)*scale);
	}

	input.update_normals();
	return input;
}
MyMesh moveVertexAloneNormal(MyMesh input, float moveDist){
	for (size_t i = 0; i < input.n_vertices(); i++)
	{
		OpenMesh::VertexHandle vi(i);
		OpenMesh::Vec3f vNorm = input.normal(vi).normalized();
		OpenMesh::Vec3f vPos = input.point(vi);
		input.set_point(vi,vPos+vNorm*moveDist);
	}
	input.update_normals();
	return input;
}
const char* filename = "../Resources/drop.obj";//bunny.obj";
const char* testModel = "D:\\WireArtExp_vc12_x64\\data\\_rabbit\\rabbit_4000f.off";
class showEdge :public script
{
private:
	glm::vec2 last_mousePos;
public:
	GLuint VBO[2];


	void onMouseEvent(int button, int state, int x, int y){
		float viewdist = glm::vec3(viewX, viewY, viewZ).length();
	}
	void onInit(){
		cout << "showModel onInit" << endl;
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		
		

		if (!OpenMesh::IO::read_mesh(mesh, meshPath, opt))
			{
				cerr << "Error: Cannot read mesh:"<<meshPath;
				return;
			}
			else
			{
				cout << "load mesh成功" << endl;
			}
		


		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}
		mesh.update_normals();


		mesh.request_vertex_texcoords2D();
		
		//cout << "点一,归一化之前为:" << mesh.point(MyMesh::VertexHandle(0));
		mesh= normalizeMesh(mesh);
		//cout << "点一,归一化之后为:" << mesh.point(MyMesh::VertexHandle(0));
		//vector<GLfloat> Vertices = vector<GLfloat>(mesh.n_faces()*9);

		Vector3f *Vertices = new Vector3f[mesh.n_faces() * 3];
		Vector3f *Normals = new Vector3f[mesh.n_faces() * 3];
		for (int f = 0; f < mesh.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				Vertices[f * 3 + vindex] = Vector3f(pos[0], pos[1], pos[2]);
				MyMesh::Normal norm = mesh.normal(*fv_it);
				Normals[f * 3 + vindex] = Vector3f(norm[0], norm[1], norm[2]);

				vindex++;
			}
			//cout << endl;
		}


		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Normals[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);


		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		//glEnable(GL_DEPTH_TEST);
		shaderHandle = CompileShaders("viewSpace.vs", "dotMap.fs");
		glUseProgram(shaderHandle);

		Vector3f vertex[6];
		vertex[0] = Vector3f(-1.0, -1.0, .0);
		vertex[1] = Vector3f(-1.0, 1.0, .0);
		vertex[2] = Vector3f(1.0, -1.0, .0);
		vertex[3] = Vector3f(-1.0, 1.0, .0);
		vertex[4] = Vector3f(1.0, -1.0, .0);
		vertex[5] = Vector3f(1.0, 1.0, .0);
		/*
		Vector2f coord[6];
		coord[0] = Vector2f(.0, 1.0);
		coord[1] = Vector2f(.0, .0);
		coord[2] = Vector2f(1.0, 1.0);
		coord[3] = Vector2f(.0, .0);
		coord[4] = Vector2f(1.0, 1.0);
		coord[5] = Vector2f(1.0, .0);*/
		Vector2f coord[6];
		coord[0] = Vector2f(.0, .0);
		coord[1] = Vector2f(.0, 1.0);
		coord[2] = Vector2f(1.0, .0);
		coord[3] = Vector2f(.0, 1.0);
		coord[4] = Vector2f(1.0, .0);
		coord[5] = Vector2f(1.0, 1.0);
		//第二个pass的资料设置
		GLuint VBO_pass2[2];
		glGenVertexArrays(1, &VAO_pass2);
		glBindVertexArray(VAO_pass2);

		glGenBuffers(2, VBO_pass2);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, 6 * sizeof(Vector3f), &vertex[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[1]);
		glBufferData(GL_ARRAY_BUFFER, 6 * sizeof(Vector2f), &coord[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_pass2[1]);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);


		//颜色buffer

		glGenTextures(1, &textureID);

		// "Bind" the newly created texture : all future texture functions will modify this texture
		//glActiveTexture(GL_TEXTURE0+100);
		glBindTexture(GL_TEXTURE_2D, textureID);
		//glActiveTexture(GL_TEXTURE0 + 100);
		// Give the image to OpenGL
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_FLOAT, NULL);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		//创建framebuffer

		glGenFramebuffers(1, &fbHandle);
		glBindFramebuffer(GL_FRAMEBUFFER, fbHandle);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);
		/*
		glGenTextures(1, &zTextureID);
		glBindTexture(GL_TEXTURE_2D, zTextureID);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_FLOAT, NULL);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, zTextureID, 0);
		*/

		//GLuint depthTextureID;
		glGenTextures(1, &depthID);
		glBindTexture(GL_TEXTURE_2D, depthID);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		//glTexStorage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, SCR_WIDTH, SCR_HEIGHT);
		//glTexStorage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32F, SCR_WIDTH, SCR_HEIGHT);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, SCR_WIDTH, SCR_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthID, 0);

		//GLenum bufs[] = {GL_COLOR_ATTACHMENT0};
		//	glDrawBuffers(1, bufs);

		cerr << glCheckFramebufferStatus(GL_FRAMEBUFFER);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		{
			cout << "frame buffer grejkgojeojop!" << endl;
		}


		glBindTexture(GL_TEXTURE_2D, textureID);

		cout << "開始編譯pass2 shader"<<endl;
		pass2_shader = CompileShaders("tex.vs", "vp_sobel.fs");

		//glUseProgram(pass2_shader);
		cout << "准备完成" << endl;
	}
	void onRender(){
		//cout << "render...";
		float speed = 0.01;

		//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		glBindFramebuffer(GL_FRAMEBUFFER, fbHandle);
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glBindVertexArray(VAO);
		glUseProgram(shaderHandle);
		//glUseProgram(shaderHandle);
		angle += speed;
		//offset -= speed;
		//計算相機位置---------------------------------------------------------------
		glm::vec4 campos(viewX,viewY,viewZ,0);
		glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f),angle,glm::vec3(0,1,0));
		glm::vec4 finpos4 = campos*rotat_y;
		glm::vec3 finpos = glm::vec3(finpos4);

		if (autoRotat){
			glm::vec4 autoPos = glm::vec4(0, 0, 2, 0)*rotat_y;
			viewX = autoPos.x;
			viewY = autoPos.y;
			viewZ = autoPos.z;
		}
		//---------------------------------------------------------------------------
		glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0, 1, 0));
		glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, offset));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 0, 0, z }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		glm::mat4 fakeOrthoMat = glm::perspective(1.0f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)*glm::lookAt(glm::vec3{ 0, 0, z }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		glm::mat4 dynFakeOrthoMat = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)*glm::lookAt(glm::vec3(viewX, viewY, viewZ), glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });

		float aspect = SCR_WIDTH / SCR_HEIGHT;
		glm::mat4 ortho = glm::ortho(-1 * aspect, 1 * aspect, -1.0f, 1.0f, nearZ, farZ);
		/*
		cout << "正交投影矩陣:" << endl;
		for (int y = 0; y < 4; y++){
			for (int x = 0; x < 4; x++){
				cout << ortho[y][x] << ", ";
			}
			cout << endl;
		}*/
		glm::mat4 orthoMat = glm::ortho(-1 * aspect, 1 * aspect, -1.0f, 1.0f, nearZ, farZ)* glm::lookAt(glm::vec3{ viewX, viewY, viewZ }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		glm::mat4 MVMat = glm::lookAt(glm::vec3{ viewX, viewY, viewZ }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		GLuint loc = glGetUniformLocation(shaderHandle, "gWVP");
		assert(loc != 0xFFFFFFFF);
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(dynFakeOrthoMat*glm::mat4(1.0f)));
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(orthoMat*glm::mat4(1.0f)));
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(fakeOrthoMat*translate*rotate));//*rotate));

		GLuint loc_MV = glGetUniformLocation(shaderHandle, "gMV");
		glUniformMatrix4fv(loc_MV,1,GL_FALSE,glm::value_ptr(MVMat*glm::mat4(1.0f)));
		
		glm::vec3 lightDir = glm::vec3(1, 0, 0);
		GLuint loc_light = glGetUniformLocation(shaderHandle, "light_direction");
		//assert(loc_light != 0xFFFFFFFF);
		glUniform3f(loc_light, lightDir[0], lightDir[1], lightDir[2]);
		//glUniform4f()
		//glClear(GL_COLOR_BUFFER_BIT);

		glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces() * 9);
		glBindVertexArray(0);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
		glDisable(GL_DEPTH_TEST);
		glClear(GL_COLOR_BUFFER_BIT);
		glUseProgram(pass2_shader);
		glBindVertexArray(VAO_pass2);


		//glBindFramebuffer(GL_FRAMEBUFFER,fbHandle);
		//glBindTexture(GL_TEXTURE_2D, textureID);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, textureID);
		GLuint loc_tex = glGetUniformLocation(pass2_shader, "ourTexture");
		glUniform1i(loc_tex, 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, depthID);
		GLuint loc_tex2 = glGetUniformLocation(pass2_shader, "depthTexture");
		glUniform1i(loc_tex2, 1);
		//assert(loc_tex != 0xFFFFFFFF);

		GLuint loc_near = glGetUniformLocation(pass2_shader, "zNear");
		glUniform1f(loc_near, nearZ);
		GLuint loc_far = glGetUniformLocation(pass2_shader, "zFar");
		glUniform1f(loc_far, farZ);

		GLuint loc_viewPos = glGetUniformLocation(pass2_shader, "viewPos");
		glm::vec3 toCam = glm::vec3(dynFakeOrthoMat*glm::vec4(viewX, viewY, viewZ,1));
		glUniform3f(loc_viewPos, toCam.x,toCam.y,toCam.z);

		GLuint loc_th = glGetUniformLocation(pass2_shader, "threshold");
		glUniform1f(loc_th, threshold);
		//cout << "toCam:" << toCam.x << "," << toCam.y << "," << toCam.z << endl;
		glDrawArrays(GL_TRIANGLES, 0, 6);
		//glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		glBindVertexArray(0);
		//z += speed;
		//glEnableVertexAttribArray(0);
		//glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
	void updateThreshold(float offset)
	{
		threshold += offset;
	}
	showEdge(){
		meshPath = "D:\\WireArtExp_vc12_x64\\data\\_rabbit\\rabbit_4000f.off";
		autoRotat = true;
	}
	showEdge(char* path){
		cout << "設置path:" << path<<endl;
		meshPath = path;
		autoRotat = true;
	}
	showEdge(float x,float y,float z){
		viewX = x;
		viewY = y;
		viewZ = z;
		meshPath = "../Resources/bunny.obj";
	}
	showEdge(float x, float y, float z, char* path){
		viewX = x;
		viewY = y;
		viewZ = z;
		meshPath = path;
		cout << "初始化showEdge meshPath:" << meshPath << endl;
	}
	/*showEdge(float x, float y, float z, char* path, char* save_path):showEdge(x,y,z,path){
		saveImgPath = save_path;
		cout << "初始化showEdge save_path:" << save_path << endl;
	}*/

private:
	MyMesh mesh;
	GLuint shaderHandle;
	GLuint pass2_shader;
	GLuint frameTexHandle;
	GLuint frameDepthHandle;
	float z = 1;
	float offset = -1;
	float angle = 0;
	float nearZ = 1.0f;
	float farZ = 60.0f;
	GLuint VAO;
	GLuint VAO_pass2;
	GLuint fbHandle;
	GLuint textureID;
	GLuint depthID;
	GLuint zTextureID;
	float viewX = 0;
	float viewY = 0;
	float viewZ = 2;
	float threshold = 0.35;
	char* meshPath;
	char* saveImgPath;
	bool autoRotat = false;
};
Vector3f toVector3f(OpenMesh::Vec3f p){
	return Vector3f(p[0], p[1], p[2]);
}
Vector3f* toVector3Arrays(vector<PointCurve> lines,vector<int> &sizes){//把lines轉化為邊端點的集合,sizes表明每一組邊的點集合有幾個點
	vector<Vector3f> *points=new vector<Vector3f>();
	for (int i = 0; i < lines.size(); i++){
		//cout << "i=" << i << ":";
		//cout << "points i=" << i << " size:" << lines[i].points.size()<<" ";
		for (int j = 0; j < lines[i].points.size()-1; j++){
			OpenMesh::Vec3f point1= lines[i].points[j].pos;
			Vector3f p1 = Vector3f(point1[0], point1[1], point1[2]);
			//cout << "point1:" << p1.x<<","<<p1.y<<","<<p1.z<<endl;
			points->push_back(Vector3f(point1[0],point1[1],point1[2]));
			OpenMesh::Vec3f point2 = lines[i].points[j+1].pos;
			//cout << "point2:" << Vector3f(point2[0], point2[1], point2[2]);
			points->push_back(Vector3f(point2[0], point2[1], point2[2]));
			
			//points.push_back(Vector3f(point[0],point[1],point[2]));
		}
		
		OpenMesh::Vec3f point_end = lines[i].points.back().pos;
		points->push_back(Vector3f(point_end[0],point_end[1],point_end[2]));
		OpenMesh::Vec3f point_start = lines[i].points.front().pos;
		points->push_back(Vector3f(point_start[0], point_start[1], point_start[2]));
		cout << endl;
		sizes.push_back(2*(lines[i].points.size()));
	}
	return points->data();
}
vector<Vector3f> indexArrayToVector3s(vector<int> indexArray,MyMesh mesh){
	vector<Vector3f> Vec3fArray = vector<Vector3f>(indexArray.size());
	for (int i = 0; i < indexArray.size(); i++){
		int index = indexArray[i];
		float* pos= mesh.point(OpenMesh::VertexHandle(index)).data();
		//cout << "pos:(" << pos[0] << "," << pos[1] << "," << pos[2] << ")->";
		Vec3fArray[i] = Vector3f(pos[0],pos[1],pos[2]);
		//cout << "Vec3fArray[" << i << "]:(" << Vec3fArray[i][0] << "," << Vec3fArray[i][1] << "," << Vec3fArray[i][2] << ")"<<endl;
		cout << endl;
	}
	return Vec3fArray;
}
vector<Vector3f> indexArrayToVector3s(vector<int> indexArray, omMesh mesh){
	vector<Vector3f> Vec3fArray = vector<Vector3f>(indexArray.size());
	for (int i = 0; i < indexArray.size(); i++){
		int index = indexArray[i];
		DGPC::Vector3<double> pos = mesh.point(OpenMesh::VertexHandle(index));
		//cout << "pos:(" << pos[0] << "," << pos[1] << "," << pos[2] << ")->";
		Vec3fArray[i] = Vector3f(pos[0], pos[1], pos[2]);
		//cout << "Vec3fArray[" << i << "]:(" << Vec3fArray[i][0] << "," << Vec3fArray[i][1] << "," << Vec3fArray[i][2] << ")"<<endl;
		cout << endl;
	}
	return Vec3fArray;
}
static bool show_demo_window = true;
string sphere_path = "../Resources/bunny.obj"; //"../Resources/uvSphere_16.obj"; //"../Resources/fork.off"; //"../Resources/raptor_2500f.obj";// "../Resources/rabbit_4000f_remesjed.obj";//"../Resources/drop.obj";//"../Resources/dog.obj";//"../Resources/hand.obj";//"../Resources/horsehead_smooth.obj"; 
//"../Resources/cube30.obj";//"../Resources/cylinder.obj";
//"../Resources/hand_remeshed.obj";//"../Resources/dog.off";//"../Resources/rabbit_4000f_remesjed.obj"; //"../Resources/uv3_fix.obj";
class showPolar :public script
{
protected:
	GLuint VAO[1];
	GLuint e_VAO;
	GLuint v_VAO;
	GLuint c_VAO;
	GLuint axis_VAO;
	//GLuint subf_VAO;
	int c_totalNum;
	polarMap* pMap;
	//polarMap_angleProjection *pMap;
	int center = 296;//406;//315;//296;//219;//2950;//1030; //82;//407;//714;//385;//684;//158;//209;
	//2950是dog的起點//1030是hand的起點//407是shark的起點 //296 是bunny的起點//315是rabbit的起點//385是cylinder的起點//315是horsehead_smooth的起點
	//406是raptor的起點
	int contourNum = 20;//10;//30;//50;
	int elemPerRing = 30;//30;
	int sampleRate = 10;
	//const int AXIS_VERTEX_NUM = 17;
	//int* zeroAxis = new int[AXIS_VERTEX_NUM]{296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 205};
	int elemId = 0;
protected:
	float control_rotate_y = 0;
	float control_rotate_x = 0;
	float rotateSpeed = 0.05f;//0.0002f;
	float control_cameraZ = 2.0f;
	float cameraSpeed = 0.01f;
public:
	MyMesh  mesh;

	GLuint vbo[2];
	GLuint e_vbo;
	GLuint v_vbo;
	GLuint c_vbo[2];
	GLuint axis_vbo[2];
	//GLuint subf_vbo[2];

	GLuint shaderHandle;
	GLuint e_shaderHandle;
	GLuint v_shaderHandle;
	GLuint p_shaderHandle;
	GLuint axis_shaderHandle;
	float nearZ = 0.5f;
	float farZ = 60.0f;

	//float rotateAngle = 0;

	vector<Vector3f> showPoints;
	bool showFace = true;//false;
	bool showEdge = true;//false;
	bool showPoint = false;
	bool showContour = false;
	bool showAxis = false;

	bool contourOnly = false;
	bool drawGDPCRespond = true;
	bool drawNodeGraph = true;
	bool drawCompositeResult = true;
	bool createMesh = true;
	bool gridOnly = false;
	bool drawZeroAxis = false;
	bool drawDivergenceConnect=false;

	bool demoSequences = false;
	int nowStep=0;
	vector<Vector3f> showEdges;
	vector<Vector4f> axisColor;
	vector<Vector4f> sequenceColors;

	vector<Vector3f> showGrids;
	vector<Vector4f> gridColors;
	vector<Vector3f> zeroAxis;
	vector<Vector4f> zeroAxisColors;
	vector<Vector3f> divergenceConnects;
	vector<Vector4f> connectColors;

	vector<Vector4f> contourColor;
	vector<Vector4f> contourColor_angle;
	float *whitePercentage;
	const int DEMO_POINT_NUM_PER_CONTOLUR=16;

	void onInit(){
		
		//使用者輸入-----------------------------------------------------------------------------
		string img_path = "../Resources/elem8.png";
		cout << "請輸入圖元的路徑:" << endl;
		cin >> img_path;
		cout << "請輸入圖元採樣頻率:" << endl;
		string sampleRate_str;
		cin >> sampleRate_str;
		sampleRate = stoi(sampleRate_str);
		cout << "請輸入3D模型的路徑:" << endl;
		cin >> sphere_path;
		cout << "請輸入起始點索引值:" << endl;
		string center_str;
		cin >> center_str;
		center = stoi(center_str);
		cout << "請輸入等高線環數:" << endl;
		string contourNum_str;
		cin >> contourNum_str;
		contourNum = stoi(contourNum_str);
		cout << "請輸入每層等高線上的圖元數:" << endl;
		string elemPerRing_str;
		cin >> elemPerRing_str;
		elemPerRing = stoi(elemPerRing_str);
		clock_t method_start, method_end;
		method_start = clock();
		//圖元處理--------------------------------------------------------------------------------
		imageConvert convert(img_path);//使用路徑讀取圖元並且創建pixelGraph
		int totalPixelNum = convert.totalPixelPointNum;//picelGraph節點的數量,並不是圖片的width*height
		vector<vector<pixelPoint>> elemLines = convert.sample(sampleRate);//採樣圖片,將graph變為以線為主,vector<pixelPoint>為一條pixelPoint組成的線,elemLines為線的陣列
		/*cout << "總共有" << elemLines.size() << "條elemLine:";
		for (size_t i = 0; i < elemLines.size(); i++)
		{
			cout << "line:" << i << ":";
			for each (pixelPoint point in elemLines[i])
			{
				cout << "(" << point.x << "," << point.y << "),";
			}
			cout << endl;
		}*/
		convert.showGraphicConnect(elemLines);//橫向匹配縱向匹配以及使用mtsp規劃圖元的連接順序
		//讀取模型------------------------------------------------------------------------------
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();

		
		if (!OpenMesh::IO::read_mesh(mesh, sphere_path, opt))//OpenMesh讀取模型
		{
			cerr << "Error: Cannot read mesh:" << sphere_path;
			return;
		}
		else
		{
			cout << "load mesh成功" << endl;
		}
		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}
		mesh.update_normals();
		mesh.request_vertex_texcoords2D();
		mesh = normalizeMesh(mesh);//歸一化模型,首先算出模型中心點,然後使離中心點最遠的頂點距離為1
		MyMesh mesh4Rending = moveVertexAloneNormal(mesh, -0.0075f);//實際渲染的mesh數據,這麼做的原因是除了直接渲染出表面之外,我們還要用mesh for rending做內襯渲染sequence和path graph,mesh4Rending會比實際上用來計算graph的mesh小一點點
		//創建pMap和contours-----------------------------------------------------------------------------
		pMap =new polarMap_angleProjection(mesh, center);//創建geodesic polar map
		pMap->meshPath = sphere_path;
		pMap->createMap();

		polarPoint farthest = pMap->getPointFrom(center);
		for (int i = 1; i < pMap->pointNum(); i++){
			if (farthest.distance < pMap->getPointFrom(i).distance){//如果第i個點距離比farthest遠則更新farthest
				farthest = pMap->getPointFrom(i);
			}
		}
		vector<Vector3f> Vertices(mesh.n_faces() * 3);
		vector<Vector3f> edge_Vertex(mesh.n_faces()*6);
		float *distPercentage = new float[mesh.n_faces() * 3];
		whitePercentage = new float[mesh.n_faces() * 3];
		vector<PointCurve> lines = pMap->getContourLine(contourNum);//創建等高線

		vector<vector<int>> axises = pMap->calBaseAxis();//計算角度
		cout << "印出contour1 angles:";
		for each (pointFromTwoSource pt in pMap->contours->operator[](1)->pfts)
		{
			cout << pt.angleByCal << ",";
		}
		cout << endl;
		if (!contourOnly)
			pMap->assignDivergence();//分配分歧等高線之間的表面
		//將axis的資料存起來之後給openGL繪製-----------------------------------------------------------
		for each (vector<int> axis in axises)
		{
			for (int v = 0; v < axis.size() - 1; v++){
				OpenMesh::Vec3f p1 = mesh.point(OpenMesh::VertexHandle(axis[v]));
				zeroAxis.push_back(Vector3f(p1[0], p1[1], p1[2]));
				zeroAxisColors.push_back(Vector4f(1.0, 0.2, 0.2, 1.0));
				OpenMesh::Vec3f p2 = mesh.point(OpenMesh::VertexHandle(axis[v+1]));
				zeroAxis.push_back(Vector3f(p2[0], p2[1], p2[2]));
				zeroAxisColors.push_back(Vector4f(1.0, 0.2, 0.2, 1.0));
			}
		}
		//將mesh的表面和邊存起來之後供openGL繪製-------------------------------------------------------------
		Vector4f color_half_black(0.3,0.3,0.3,1);
		Vector4f color_default(0.2, 0.2, 0.8, 1);
		Vector4f colorTable[4]{Vector4f(1.0,0.2,0.2,1.0),Vector4f(0.2,1.0,0.2,1.0),Vector4f(0.2,0.2,1.0,1.0),Vector4f(0.0,0.0,0.0,1.0)};

		for (int f = 0; f < mesh4Rending.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;

			for (MyMesh::FaceVertexIter fv_it = mesh4Rending.fv_begin(fh); fv_it != mesh4Rending.fv_end(fh); fv_it++)
			{
				float* pos = mesh4Rending.point(*fv_it).data();
				Vertices[f * 3 + vindex] = Vector3f(pos[0], pos[1], pos[2]);
				//distPercentage[f * 3 + vindex] = pMap->getPointFrom((*fv_it).idx()).distance/farthest.distance;
				distPercentage[f * 3 + vindex] = pMap->getPercentageFromVertex((*fv_it).idx()); //pMap->getPointFrom((*fv_it).idx()).angle / (2 * M_PI);
				whitePercentage[f * 3 + vindex] = 0.9f;
				vindex++;
			}
			if (f < 10){
				cout << endl;
			}
			//第二個pass的邊
			int eindex = 0;
			for (MyMesh::FaceHalfedgeCWIter fhe_cit = mesh4Rending.fh_cwbegin(fh); fhe_cit != mesh4Rending.fh_cwend(fh); fhe_cit++){
				float* from_pos = mesh4Rending.point(mesh4Rending.from_vertex_handle(*fhe_cit)).data();
				float* to_pos = mesh4Rending.point(mesh4Rending.to_vertex_handle(*fhe_cit)).data();
				//cout << "index:" << f * 6 + eindex * 2 << "from pos:(" << from_pos[0] << "," << from_pos[1] << "," << from_pos[2] << ")";
				edge_Vertex[f * 6 + eindex * 2] = Vector3f(from_pos[0], from_pos[1], from_pos[2]);
				//cout << "index" << f * 6 + eindex * 2 +1<< "to pos:(" << to_pos[0] << "," << to_pos[1] << "," << to_pos[2] << ")" << endl;
				edge_Vertex[f * 6 + eindex * 2 + 1] = Vector3f(to_pos[0], to_pos[1], to_pos[2]);
				eindex++;
			}
			//cout << endl;
		}
		showPoints.clear();
		MyMesh::VertexHandle startPoint = MyMesh::VertexHandle(center);
		OpenMesh::Vec3f point = mesh.point(startPoint);
		showPoints.push_back(Vector3f(point[0], point[1], point[2]));

		const int debugElemIdx = 15;
		const int debugELineIdx = 0;
		//計算表面的等分-------------------------------------------------------------------------------
		bool testGDCPath = false;
		optimizerGroup optGroup(*(pMap->contours), elemPerRing);
		optGroup.stage_solve(pMap->contours->operator[](0), convert.Shift_percent);//等分表面
		/*cout << "stage solve完成."<<endl;
		for (size_t i = 0; i < pMap->contours->size(); i++)
		{
			offsetInf inf = optGroup.origenInf[i];
			cout << "第" << i << "offsetInf: origen contour:" << inf.origenContourIndex << "level:" << inf.subLevel << "total_offset:" << inf.totalOffset << ";"<<endl;
			for each (float angle in inf.origenAngleList)
			{
				cout << angle << ",";
			}
			cout << endl;
		}
		float angle59= optGroup.getXat(5, 9, 0.5);
		cout << "optGroup.getXat(5, 9, 0.5):" << angle59 << endl;*/
		//放置圖元到等高線之間(非分歧的部分)----------------------------------------------------------------
		int debugContourIndexList[1] = { 13 };//{ 5, 8, 11, 15, 17, 18, 19, 20, 21, 22 };
		GraphArray nGraph(pMap->contours->size());//這是pathGraph,之後會不斷添加node到這個graph里
		const float NGROUP_Y_FIX = 0.99f;//當percent_y = 1.0的時候 percent_y+i會跳到contour i+1的percent_y=0去,但是i的下一條contour 不一定是i+1
		//所以需要乘一個小於1的NGROUP_Y_FIX來區分i的1.0 和i+1的0.0
		vector<int> contour_failCounter(pMap->contours->size());
		//for (int i = 0; i < 0; i++)
		if (!contourOnly){
			clock_t total_start_time,total_end_time;
			total_start_time = clock();
			showEdges.clear();
			axisColor.clear();
			cout << "創建pathGraph,非分歧的部分."<<endl;
			
			for (int i = 0; i < pMap->contours->size(); i++)
			{
				contour_failCounter[i] = 0;
				if (pMap->contours->operator[](i)->nextContours.size() != 1){
					continue;
				}
				cout << "contour:" << i <<":";
				clock_t start_t,end_t;
				start_t = clock();

				vector<float> last_x;
				vector<float> next_x;
				int i_next = pMap->contours->operator[](i)->nextContours[0]->id;

				last_x = optGroup.getAngleListForContour(i, true);
				last_x.push_back(last_x[0] + 1);
				next_x = optGroup.getAngleListForContour(i, false);//(i_next);
				next_x.push_back(next_x[0] + 1);

				//cout << "contour:" << pMap->contours->operator[](i)->id - 1 << ":";
				//for (int j = debugElemIdx; j < debugElemIdx + 1;j++){
				vector<float> ring1 = optGroup.getAngleListForContour(i, true);
				vector<float> ring2 = optGroup.getAngleListForContour(i, false);
				int failNum = 0;
				for (int j = 0; j < ring1.size(); j++){//第j個間隔

					//cout << "畫出 debugElemIdx:" << j << " elemId:" << elemId << endl;
					if (i==5)
						cout << "j=" << j;
					if (true){
						contour *nowC = pMap->contours->operator[](i);
						contour *nextC = nowC->nextContours[0];
						float x1 = optGroup.getXat(i, j, 0);//ring1[j];
						float x2 = optGroup.getXat(i, j, 0);//ring2[j];
						//cout << x2 << ", ";
						float c = (float)j / (float)elemPerRing;
						OpenMesh::Vec3f p1 = nowC->getPosAtAngle(x1);
						showGrids.push_back(Vector3f(p1[0], p1[1], p1[2]));
						//axisColor.push_back(Vector4f(0.3,1,0.3,1));
						//axisColor.push_back(Vector4f(c, c, c, 1));
						gridColors.push_back(colorTable[2]);
						OpenMesh::Vec3f p2 = nextC->getPosAtAngle(x2);
						showGrids.push_back(Vector3f(p2[0], p2[1], p2[2]));
						//axisColor.push_back(Vector4f(0.3, 1, 0.3, 1));
						//axisColor.push_back(Vector4f(c, c, c, 1));
						gridColors.push_back(colorTable[2]);
					}
					if(true){
						bool elemFailed = false;
						//vector<Vector3f> tempPoint;
						//vector<Vector4f> tempColor;
						int next_i = pMap->contours->operator[](i)->nextContours[0]->id;
						//float maxAngle1 = ring1.back();
						//float maxAngle2 = ring2.back();
						vector<path> tempPaths;//臨時記錄圖元需要繪製的路徑,如果在繪製投影路徑中有一次失敗失敗,則清空tempPaths,重新使用測地線繪製整個圖元

						for each (vector<pixelPoint> eline in elemLines)
						{
							path epath;//用於組GraphArray
							//MyMesh::FaceHandle f1 = pMap->contours->operator[](i)->getFaceAtAngle(mesh);
							//OpenMesh::Vec3f n1= mesh.normal(f1);

							if (eline.front().percentage_y != 1.0f){//初始化epath
								/*epath.start_2d[0] = last_x[j] + (last_x[j + 1] - last_x[j])*eline.front().percentage_x;
								if (epath.start_2d[0] > maxAngle1){//為什麼要有這個修正的原因是last_x在產生的時候把最前端的角度拼接到了尾端來產生一個遞增結果
									//這總結果在繪製的時候好用,但是在產生ngraph的時候就會出現最尾端的點被解讀成了兩種角度,所以要把所有大於原最尾端的點修正會第一點
									epath.start_2d[0] -= 1;
								}*/
								epath.start_2d[0] = optGroup.getXat(i,j,eline.front().percentage_x);//j + eline.front().percentage_x;
								epath.start_2d[1] = eline.front().percentage_y*NGROUP_Y_FIX + i;
								/*if (epath.start_2d[1] == 9){
									cout << "start y==" << 9 << " next_x[j+1]:" << next_x[j + 1] << " next_x[j]:" << next_x[j] << " interval:" << next_x[j + 1] - next_x[j] << " offset x:" << (next_x[j + 1] - next_x[j])*eline.back().percentage_x<<"final x:"<<epath.start_2d[0];
								}*/
								epath.pixelId[0] = eline.front().id + elemId*totalPixelNum;//+elemId*totalPixelNum保證每個不重合的端點pixelId不重複
							}
							else{//如果percent_y ==1的話,說明這個點實際上是下一條contour的點

								/*epath.start_2d[0] = next_x[j] + (next_x[j + 1] - next_x[j])*eline.front().percentage_x;
								if (epath.start_2d[0] > maxAngle2){
									epath.start_2d[0] -= 1;
								}*/
								epath.start_2d[0] = optGroup.getXat(i, j, eline.front().percentage_x);
								epath.start_2d[1] = next_i;
								epath.pixelId[0] = eline.front().id + elemId*totalPixelNum;
							}

							if (eline.back().percentage_y != 1.0f){
								/*epath.end_2d[0] = last_x[j] + (last_x[j + 1] - last_x[j])*eline.back().percentage_x;
								if (epath.end_2d[0] > maxAngle1){
									epath.end_2d[0] -= 1;
								}*/
								epath.end_2d[0] = optGroup.getXat(i, j, eline.back().percentage_x);
								epath.end_2d[1] = eline.back().percentage_y*NGROUP_Y_FIX + i;
								epath.pixelId[1] = eline.back().id + elemId*totalPixelNum;
							}
							else{

								/*epath.end_2d[0] = next_x[j] + (next_x[j + 1] - next_x[j])*eline.back().percentage_x;
								if (epath.end_2d[0] > maxAngle2){
									epath.end_2d[0] -= 1;
								}*/
								epath.end_2d[0] = optGroup.getXat(i, j, eline.back().percentage_x);
								epath.end_2d[1] = next_i;
								if (epath.end_2d[1] == 9){
									cout << "end y==" << 9 << " next_x[j+1]:" << next_x[j + 1] << " next_x[j]:" << next_x[j] << " interval:" << next_x[j + 1] - next_x[j] << " offset x:" << (next_x[j + 1] - next_x[j])*eline.back().percentage_x << "final x:" << epath.end_2d[0];
								}
								epath.pixelId[1] = eline.back().id + elemId*totalPixelNum;
							}
							//cout << "draw epath start(" << epath.start_2d[0] << "," << epath.start_2d[1] << ") => end(" << epath.end_2d[0] << "," << epath.end_2d[1] << ")"<<endl;
							bool fail = false;
							for (int e = 0; e < eline.size() - 1; e++)//eline是在曲線按順序上取樣后的到的點集合
								//for (int e = debugELineIdx; e < debugELineIdx + 1; e++)
							{//因為畫的是線段所以每次需要兩個點:起點和終點
								//cout << "e:" << e;

								//float x1 = last_x[j] + eline[e].percentage_x*(last_x[(j + 1)] - last_x[j]);//eline[e].percentage_x*0.05 + j;
								float x1 = optGroup.getXat(i, j, eline[e].percentage_x);
								/*float x12 = 0;

								if (next_x[j] > next_x[j + 1]){//說明跨過了0度軸

									cout << "next_x[j]:" << next_x[j] << " > next_x[j + 1]:" << next_x[j + 1];
									system("pause");
								}
								else{
									x12 = next_x[j] + eline[e].percentage_x*(next_x[(j + 1)] - next_x[j]);
								}*/
								float x12 = optGroup.getXat(i, j, eline[e].percentage_x);
								float y1 = eline[e].percentage_y;
								//if (i==5)
								//	cout << "(x1:" << x1 <<",x12:"<<x12<< ",y1:" << y1<<")";
								OpenMesh::Vec3f p1 = pMap->getPointAtBtw(i, x1, x12, y1, fail);//起點由x1 x12內插出來
								epath.point3ds.push_back(p1);
								//cout << "projection p1:" << p1 << " ";
								if (fail){
									//cout << "失敗 ";
									failNum++;
									elemFailed = true;
									break;
								}
								//tempPoint.push_back(Vector3f(p1[0], p1[1], p1[2]));
								//Vector4f c1(j*0.05, j*0.05, j*0.05, 1);
								//Vector4f c1(0.25f, 0.25f, 0.25f, 1);
								//tempColor.push_back(color_half_black);
								//else{
								//float x2 = last_x[j] + eline[e + 1].percentage_x*(last_x[(j + 1)] - last_x[j]);
								/*float x22 = 0;//next_x[j] + eline[e + 1].percentage_x*(next_x[(j + 1)] - next_x[j]);
								if (next_x[j] > next_x[j + 1]){//說明跨過了0度軸
									system("pause");
								}
								else{
									x22 = next_x[j] + eline[e + 1].percentage_x*(next_x[(j + 1)] - next_x[j]);
								}*/
								float x2 = optGroup.getXat(i,j,eline[e+1].percentage_x);
								float x22 = optGroup.getXat(i, j, eline[e + 1].percentage_x);
								float y2 = eline[e + 1].percentage_y;
								//if(i==5)
								//	cout << "(x2:" << x2 << "x22:" << x22 << ",y2:" << y2 << ")";
								OpenMesh::Vec3f p2;

								//total_contour_num++;

								p2 = pMap->getPointAtBtw(i, x2, x22, y2, fail);
								if (e == eline.size() - 2){
									epath.point3ds.push_back(p2);
								}
								//cout << "projection p2:" << p2 << "; ";
								if (fail){
									//cout << "失敗 ";
									failNum++;
									elemFailed = true;
									break;
								}


							}
							if (fail){
								break;
							}
							tempPaths.push_back(epath);
						}
						if (elemFailed){//投影路徑失敗改為申請測地線
							//tempPoint.clear();
							//tempColor.clear();
							tempPaths.clear();
							//map<float, vector<OpenMesh::Vec3f>> pathRecords;
							for each (vector<pixelPoint> eline in elemLines)
							{
								vector<pathRequestPoint> lineRquest;
								for each (pixelPoint pt in eline)
								{
									float y1 = pt.percentage_y;
									//float x1 = last_x[j] + pt.percentage_x*(last_x[(j + 1)] - last_x[j]);//eline[e].percentage_x*0.05 + j;
									/*float x12 = 0;
									if (next_x[j] > next_x[j + 1]){//說明跨過了0度軸
										float temp_x = next_x[j] - 1.0f;
										x12 = temp_x + pt.percentage_x*(next_x[(j + 1)] - temp_x);
										if (x12 < 0){
											x12 = -x12;
										}
									}
									else{
										x12 = next_x[j] + pt.percentage_x*(next_x[(j + 1)] - next_x[j]);
									}*/
									float x1 = optGroup.getXat(i, j, pt.percentage_x);
									float x12 = optGroup.getXat(next_i,j,pt.percentage_x);
									int pixelId = pt.id + elemId*totalPixelNum;
									pathRequestPoint newRquest(i, x1, i_next, x12, y1, pixelId);
									lineRquest.push_back(newRquest);
									contour_failCounter[i]++;
									//cout << "contour_failCounter[" << i << "]:" << contour_failCounter[i];
								}
								pMap->recordGdcPointLineBtw(lineRquest);
							}
						}

						//showEdges.insert(showEdges.end(), tempPoint.begin(), tempPoint.end());
						//axisColor.insert(axisColor.end(), tempColor.begin(), tempColor.end());
						for each (path p in tempPaths)
						{
							nGraph.addAPath(p);
						}
					}
					elemId++;
				}
				end_t = clock();
				double time = ((double)end_t - start_t) / (CLOCKS_PER_SEC);
				cout << "耗時:" << time << ".共" << ring1.size()<<"elem,失敗"<<failNum<<"個"<<endl;
			}
			//放置圖元到等高線之間(分歧的部分)------------------------------------------------------------------
			cout << "創建pathGraph,分歧的部分."<<endl;
			for (int i = 0; i < pMap->contours->size(); i++){//畫出分歧的連接部分
				if (pMap->contours->operator[](i)->nextContours.size()>1){
					cout << "分歧的contour:" <<i<< ":";

					clock_t start_t, end_t;
					start_t = clock();
					contour *nowContour = pMap->contours->operator[](i);
					int count = 0;
					int failNum = 0;
					for each (vector<fp_mapping> group in nowContour->divergenceSegm)
					{

						for each (fp_mapping mapping in group)
						{
							if (mapping.gridMarkers.size()>0){ //if (mapping.intersectionMapping.size() > 1){
								
								for (int j = 0; j < mapping.gridMarkers.size();j++){//for (int j = 0; j < mapping.intersectionMapping.size() - 1; j++){
									divergenceGridInf grid = mapping.gridMarkers[j];
									if (true){
										contour *nowC = pMap->contours->operator[](i);
										contour *nextC = nowC->nextContours[count];
										//float x1 = mapping.intersectionMapping[j][0];
										float x1 = optGroup.getXat(i,grid.gridIndex,0);
										float x2 = optGroup.getXat(i, grid.gridIndex,0);//mapping.intersectionMapping[j][1];
										//cout << x2 << ", ";
										//float c = (float)j / (float)elemPerRing;
										OpenMesh::Vec3f p1 = nowC->getPosAtAngle(x1);
										showGrids.push_back(Vector3f(p1[0], p1[1], p1[2]));
										//axisColor.push_back(Vector4f(0.3,1,0.3,1));
										gridColors.push_back(colorTable[2]);
										OpenMesh::Vec3f p2 = nextC->getPosAtAngle(x2);
										showGrids.push_back(Vector3f(p2[0], p2[1], p2[2]));
										//axisColor.push_back(Vector4f(0.3, 1, 0.3, 1));
										gridColors.push_back(colorTable[2]);
									}
									//else{
									if(true){
										bool elemFailed = false;
										//vector<Vector3f> tempPoint;
										//vector<Vector4f> tempColor;
										vector<path> tempPaths;
										for each (vector<pixelPoint> eline in elemLines)
										{
											path epath;//用於組GraphArray
											if (eline.front().percentage_y != 1.0f)
											{
												epath.start_2d[0] = optGroup.getXat(i, grid.gridIndex, eline.front().percentage_x);//mapping.intersectionMapping[j][0] + eline.front().percentage_x*(mapping.intersectionMapping[(j + 1)][0] - mapping.intersectionMapping[j][0]);//last_x[j] + (last_x[j + 1] - last_x[j])*eline.front().percentage_x;
												epath.start_2d[1] = eline.front().percentage_y*NGROUP_Y_FIX + i;
												epath.pixelId[0] = elemId*totalPixelNum + eline.front().id;
											}
											else
											{
												float x12;
												//if (mapping.intersectionMapping[j][1] > mapping.intersectionMapping[j + 1][1]){//說明跨過了0度軸
												if (grid.angleCorresponds[0] > grid.angleCorresponds[1]){
													float temp_x = grid.angleCorresponds[0] - 1.0f;//mapping.intersectionMapping[j][1] - 1.0f;
													x12 = temp_x + eline.front().percentage_x*(grid.angleCorresponds[1] - temp_x);//(mapping.intersectionMapping[(j + 1)][1] - temp_x);
													if (x12 < 0){
														x12 = -x12;
													}
												}
												else{
													//x12 = mapping.intersectionMapping[j][1] + eline.front().percentage_x*(mapping.intersectionMapping[(j + 1)][1] - mapping.intersectionMapping[j][1]);
													x12 = grid.angleCorresponds[0] + eline.front().percentage_x*(grid.angleCorresponds[1]-grid.angleCorresponds[0]);
												}
												epath.start_2d[0] = x12;
												epath.start_2d[1] = pMap->contours->operator[](i)->nextContours[count]->id;
												epath.pixelId[0] = elemId*totalPixelNum + eline.front().id;
											}
											if (eline.back().percentage_y != 1.0f)
											{
												epath.end_2d[0] = optGroup.getXat(i, grid.gridIndex, eline.back().percentage_x);//mapping.intersectionMapping[j][0] + eline.back().percentage_x*(mapping.intersectionMapping[(j + 1)][0] - mapping.intersectionMapping[j][0]);//last_x[j] + (last_x[j + 1] - last_x[j])*eline.front().percentage_x;
												epath.end_2d[1] = eline.back().percentage_y*NGROUP_Y_FIX + i;
												epath.pixelId[1] = elemId*totalPixelNum + eline.back().id;
											}
											else
											{
												float x12;
												//if (mapping.intersectionMapping[j][1] >  mapping.intersectionMapping[j + 1][1]){//說明跨過了0度軸
												if (grid.angleCorresponds[0]>grid.angleCorresponds[1]){
													float temp_x = grid.angleCorresponds[0] - 1;///mapping.intersectionMapping[j][1] - 1.0f;
													x12 = temp_x + eline.back().percentage_x*(grid.angleCorresponds[1] - temp_x);//(mapping.intersectionMapping[(j + 1)][1] - temp_x);
													if (x12 < 0){
														x12 = -x12;
													}
												}
												else{
													//x12 = mapping.intersectionMapping[j][1] + eline.back().percentage_x*(mapping.intersectionMapping[(j + 1)][1] - mapping.intersectionMapping[j][1]);
													x12 = grid.angleCorresponds[0] + eline.back().percentage_x*(grid.angleCorresponds[1] - grid.angleCorresponds[0]);
												}
												epath.end_2d[0] = x12;
												epath.end_2d[1] = pMap->contours->operator[](i)->nextContours[count]->id;
												epath.pixelId[1] = elemId*totalPixelNum + eline.back().id;
											}
											bool fail = false;
											for (int e = 0; e < eline.size() - 1; e++)//eline是在曲線按順序上取樣后的到的點集合
												//for (int e = debugELineIdx; e < debugELineIdx + 1; e++)
											{//因為畫的是線段所以每次需要兩個點:起點和終點
												//cout << "e:" << e;

												float x1 = optGroup.getXat(i,grid.gridIndex,eline[e].percentage_x);//mapping.intersectionMapping[j][0] + eline[e].percentage_x*(mapping.intersectionMapping[(j + 1)][0] - mapping.intersectionMapping[j][0]);//eline[e].percentage_x*0.05 + j;
												float x12 = 0;

												//if (mapping.intersectionMapping[j][1] >  mapping.intersectionMapping[j + 1][1]){//說明跨過了0度軸
												if(grid.angleCorresponds[0]>grid.angleCorresponds[1]){
													float temp_x = grid.angleCorresponds[0] - 1.0f;//mapping.intersectionMapping[j][1] - 1.0f;
													x12 = temp_x + eline[e].percentage_x*(grid.angleCorresponds[1] - temp_x);//(mapping.intersectionMapping[(j + 1)][1] - temp_x);
													if (x12 < 0){
														x12 = -x12;
													}
												}
												else{
													x12 = grid.angleCorresponds[0] + eline[e].percentage_x*(grid.angleCorresponds[1] - grid.angleCorresponds[0]);//mapping.intersectionMapping[j][1] + eline[e].percentage_x*(mapping.intersectionMapping[(j + 1)][1] - mapping.intersectionMapping[j][1]);
												}
												float y1 = eline[e].percentage_y;
												//cout << "(x1:" << x1 <<",x12:"<<x12<< ",y1:" << y1<<")";
												OpenMesh::Vec3f p1 = pMap->getPointAtBtw(i, count, x1, x12, y1, fail);//起點由x1 x12內插出來
												epath.point3ds.push_back(p1);
												//cout << "projection p1:" << p1 << " ";
												if (fail){
													//cout << "失敗 ";
													failNum++;
													elemFailed = true;
													break;
												}

												float x2 = optGroup.getXat(i,grid.gridIndex,eline[e].percentage_x);//mapping.intersectionMapping[j][0] + eline[e + 1].percentage_x*(mapping.intersectionMapping[(j + 1)][0] - mapping.intersectionMapping[j][0]);
												float x22 = 0;//next_x[j] + eline[e + 1].percentage_x*(next_x[(j + 1)] - next_x[j]);
												//if (mapping.intersectionMapping[j][1] > mapping.intersectionMapping[j + 1][1]){//說明跨過了0度軸
												if (grid.angleCorresponds[0]>grid.angleCorresponds[1]){
													float temp_x =grid.angleCorresponds[0] - 1.0f;
													x22 = temp_x + eline[e + 1].percentage_x*(grid.angleCorresponds[1] - temp_x);//(mapping.intersectionMapping[(j + 1)][1] - temp_x);
													if (x22 < 0){
														x22 = -x22;
													}
												}
												else{
													x22 = grid.angleCorresponds[0] + eline[e + 1].percentage_x*(grid.angleCorresponds[1] - grid.angleCorresponds[0]); //mapping.intersectionMapping[j][1] + eline[e + 1].percentage_x*(mapping.intersectionMapping[(j + 1)][1] - mapping.intersectionMapping[j][1]);
												}
												float y2 = eline[e + 1].percentage_y;

												OpenMesh::Vec3f p2;



												p2 = pMap->getPointAtBtw(i, count, x2, x22, y2, fail);
												if (e == eline.size() - 2){
													epath.point3ds.push_back(p2);
												}
												if (fail){
													//cout << "失敗 ";
													failNum++;
													elemFailed = true;
													break;
												}
											}
											if (fail){
												break;
											}
											tempPaths.push_back(epath);
										}
										if (elemFailed){//投影路徑失敗改為測地線
											//tempPoint.clear();
											//tempColor.clear();
											tempPaths.clear();
											map<float, vector<OpenMesh::Vec3f>> pathRecords;
											for each (vector<pixelPoint> eline in elemLines)
											{
												vector<pathRequestPoint> lineRquest;
												for each (pixelPoint pt in eline)
												{
													float y1 = pt.percentage_y;
													float x1 = optGroup.getXat(i, grid.gridIndex, pt.percentage_x);//mapping.intersectionMapping[j][0] + pt.percentage_x*(mapping.intersectionMapping[(j + 1)][0] - mapping.intersectionMapping[j][0]);//eline[e].percentage_x*0.05 + j;
													float x12 = 0;
													//if (mapping.intersectionMapping[j][1] > mapping.intersectionMapping[j + 1][1]){//說明跨過了0度軸
													if (grid.angleCorresponds[0]>grid.angleCorresponds[1]){
														float temp_x = grid.angleCorresponds[0] - 1.0f; //mapping.intersectionMapping[j][1] - 1.0f;
														x12 = temp_x + pt.percentage_x*(grid.angleCorresponds[1]-temp_x);//(mapping.intersectionMapping[(j + 1)][1] - temp_x);
														if (x12 < 0){
															x12 = -x12;
														}
													}
													else{
														x12 = grid.angleCorresponds[0] + pt.percentage_x*(grid.angleCorresponds[1]-grid.angleCorresponds[0]);//mapping.intersectionMapping[j][1] + pt.percentage_x*(mapping.intersectionMapping[(j + 1)][1] - mapping.intersectionMapping[j][1]);
													}
													int i_next = pMap->contours->operator[](i)->nextContours[count]->id;
													pathRequestPoint newRquest(i, x1, i_next, x12, y1, elemId*totalPixelNum + pt.id);
													lineRquest.push_back(newRquest);
													contour_failCounter[i]++;
													//cout << "contour_failCounter[" << i << "]:" << contour_failCounter[i];
												}
												pMap->recordGdcPointLineBtw(lineRquest);
											}
										}
										//showEdges.insert(showEdges.end(), tempPoint.begin(), tempPoint.end());
										//axisColor.insert(axisColor.end(), tempColor.begin(), tempColor.end());
										for each (path p in tempPaths)
										{
											nGraph.addAPath(p);
										}
									}
									elemId++;
								}
							}
						}
						count++;
					}
					end_t = clock();
					double time = ((double)end_t - start_t) / (CLOCKS_PER_SEC);
					cout << "耗時:" << time << ".失敗" << failNum << "次" << endl;
				}
			}
			total_end_time = clock();
			float total_time = ((double)(total_end_time - total_start_time)) / CLOCKS_PER_SEC;
			cout << "使用投影路徑耗時:" << total_time << "秒" << endl;
		}
		//繪製請求的測地線的部分----------------------------------------------------------------
		if (drawGDPCRespond){//
			clock_t total_start_time, total_end_time;
			total_start_time = clock();
			vector<vector<pathRequestPoint>> results = pMap->RespondsGdcPoints();
			cout << "開始請求測地線部分."<<endl;
			for each (vector<pathRequestPoint> line in results)
			{
				path newLine;
				if (line.front().y != 1.0f)
				{
					newLine.start_2d[0] = line.front().contour_angle[0];
					newLine.start_2d[1] = line.front().y*NGROUP_Y_FIX + line.front().contour_index[0];
					newLine.pixelId[0] = line.front().pixelId;
				}
				else
				{
					newLine.start_2d[0] = line.front().contour_angle[1];
					newLine.start_2d[1] = line.front().contour_index[1];//line.front().y + line.front().contour_index[0];
					newLine.pixelId[0] = line.front().pixelId;
				}
				if (line.back().y != 1.0f)
				{
					newLine.end_2d[0] = line.back().contour_angle[0];
					newLine.end_2d[1] = line.back().y*NGROUP_Y_FIX + line.back().contour_index[0];
					newLine.pixelId[1] = line.back().pixelId;
				}
				else{
					newLine.end_2d[0] = line.back().contour_angle[1];
					newLine.end_2d[1] = line.back().contour_index[1];
					newLine.pixelId[1] = line.back().pixelId;
				}
				for (size_t i = 0; i < line.size() - 1; i++)
				{
					OpenMesh::Vec3f p1 = line[i].pointPos;
					//showEdges.push_back(Vector3f(p1[0], p1[1], p1[2]));
					axisColor.push_back(color_half_black);
					newLine.point3ds.push_back(p1);

					OpenMesh::Vec3f p2 = line[i + 1].pointPos;
					showEdges.push_back(Vector3f(p2[0], p2[1], p2[2]));
					//axisColor.push_back(color_half_black);
					if (i == line.size() - 2){
						newLine.point3ds.push_back(p2);
					}
				}
				nGraph.addAPath(newLine);
				//cout << endl;
			}
			total_end_time = clock();
			cout << "使用geodesic總耗時:" << ((double)(total_end_time - total_start_time)) / CLOCKS_PER_SEC<<endl;
		}
		cout << "pathGraph創建完畢." << endl;
		//將graph儲存之後繪製-----------------------------------------------------------------------
		if (drawNodeGraph){
			showEdges.clear();
			axisColor.clear();

			vector<pathNode*> haveDraw;
			vector<pathNode*> candidates;
			candidates.push_back(nGraph.array[0][0]);
			int iter_time = 0;
			while (candidates.size() > 0)
			{
				//cout << "迭代" << iter_time++ << "; ";
				pathNode* nowNode = candidates[0];
				haveDraw.push_back(nowNode);
				candidates.erase(candidates.begin());
				for (size_t i = 0; i < nowNode->connects.size(); i++)//畫出到所有連接著的node的路徑
				{
					pathNode* next = nowNode->connects[i];
					//cout << "path size:" << nowNode->path.size();
					if (!contain(haveDraw, next)){//如果
						if (nowNode->path[i].size() == 0)
						{
							showEdges.push_back(toVector3f(nowNode->position));
							axisColor.push_back(colorTable[3]);
							showEdges.push_back(toVector3f(next->position));
							axisColor.push_back(colorTable[3]);
							//cout << "(" << nowNode->position << ")=>(" << next->position << ")";
						}
						else
						{

							showEdges.push_back(toVector3f(nowNode->position));
							axisColor.push_back(colorTable[3]);
							showEdges.push_back(toVector3f(nowNode->path[i].front()));
							axisColor.push_back(colorTable[3]);
							for (size_t j = 0; j < nowNode->path[i].size() - 1; j++)
							{
								OpenMesh::Vec3f p1 = nowNode->path[i][j];
								showEdges.push_back(toVector3f(p1));
								axisColor.push_back(colorTable[3]);
								OpenMesh::Vec3f p2 = nowNode->path[i][j + 1];
								showEdges.push_back(toVector3f(p2));
								axisColor.push_back(colorTable[3]);
							}
							showEdges.push_back(toVector3f(nowNode->path[i].back()));
							axisColor.push_back(colorTable[3]);
							showEdges.push_back(toVector3f(next->position));
							axisColor.push_back(colorTable[3]);
						}
						//if (next->y_2d < 2){
						if (!contain(candidates, next))
						{
							candidates.push_back(next);
						}
						//}
					}
				}
			}
			//將連接順序按不同的wire不同的顏色儲存在sequenceColors中--------------------------------------------------
			if (drawCompositeResult){
				showEdges.clear();
				axisColor.clear();
				vector<pair<pathNode*, pathNode*>> startDirects;
				for each (vector<pathNode*> level in nGraph.array)
				{
					for each (pathNode* node in level)
					{
						if (loopClampTo01(node->y_2d) == 0.0f)//交點
						{
							vector<pathNode*> upperNodes;
							vector<pathNode*> lowerNodes;
							for each (pathNode* neighber in node->connects)
							{
								int neighberId = neighber->pixelIds[0] % totalPixelNum;
								int elemNumber = neighber->pixelIds[0] - neighberId;
								int nodePixelId = -1;
								for each (int sumId in node->pixelIds)
								{
									if (sumId >= elemNumber&&sumId<elemNumber + totalPixelNum)
									{
										nodePixelId = sumId - elemNumber;
										break;
									}
								}
								bool isExtraStart = false;
								if (neighber->y_2d > node->y_2d)//neighbor在node之下,那麼pixel graph中,neighbor一定是start direct的to方向
								{
									if (convert.inExtraStartDirect(nodePixelId, neighberId)){
										startDirects.push_back(pair<pathNode*, pathNode*>(node,neighber));
										isExtraStart = true;
									}
								}
								else if (neighber->y_2d<node->y_2d)
								{
									if (convert.inExtraStartDirect(neighberId, nodePixelId)){
										startDirects.push_back(pair<pathNode*, pathNode*>(neighber, node));
										isExtraStart = true;
									}
								}
								if (!isExtraStart)
								{
									if (neighber->y_2d<node->y_2d)//上一層
									{
										upperNodes.push_back(neighber);
									}
									else if (neighber->y_2d>node->y_2d){
										lowerNodes.push_back(neighber);
									}
								}
							}
							if (upperNodes.size() >= lowerNodes.size()){
								for (size_t i = 0; i < lowerNodes.size(); i++)
								{
									node->directedConnects[upperNodes[i]] = lowerNodes[i];
								}
							}
							else{
								//cout << "node:(" << node->x_2d << "," << node->y_2d << "):upperNodes:" << upperNodes.size() << "lowerNodes:" << lowerNodes.size() << endl;
								for (size_t i = 0; i < upperNodes.size(); i++)
								{
									node->directedConnects[upperNodes[i]] = lowerNodes[i];
								}
								for (size_t i = upperNodes.size(); i < lowerNodes.size(); i++)
								{
									startDirects.push_back(pair<pathNode*, pathNode*>(node, lowerNodes[i]));
								}
							}
						}
						else{//圖元內部點
							for each (int sumId in node->pixelIds)//由於一個點可能在不同的圖元中扮演不同的像素點,pixelIds會是一個陣列
							{

								int pixelId = sumId%totalPixelNum;
								int elemNumber = sumId - pixelId;
								map<int, int> mapping = convert.connectDirect[pixelId];
								vector<pair<pathNode*, int>> pixelLevelConnect;//先找到周圍的點中有哪些是在像素圖的對應點中有相連
								for each (pathNode* neighbor in node->connects)
								{
									//cout << "測試neighbor(" << neighbor->x_2d << "," << neighbor->y_2d << "):";
									for each (int pId in neighbor->pixelIds)
									{
										//cout << "pId:" << pId;
										if (pId >= elemNumber&&pId < elemNumber + totalPixelNum){//用pId和elemNumber的關係來測試neiber有沒有像素圖內部的關係
											//cout << " connect as:" << pId - elemNumber << ";";
											pixelLevelConnect.push_back(pair<pathNode*, int>(neighbor, pId - elemNumber));
										}
									}
								}
								for each (pair<pathNode*, int> pair1 in pixelLevelConnect)
								{
									int lastId = pair1.second;
									if (convert.inExtraStartDirect(pixelId, lastId)){
										//cout << "startDirects set!";
										startDirects.push_back(pair<pathNode*, pathNode*>(node, pair1.first));
									}
									if (mapping.find(lastId) == mapping.end()){//在mapping中找不到,大概因為lastId是一個終點吧,mapping的key是起點value是終點,
										continue;
									}
									int nextId = mapping[lastId];//推算出下一個節點的pixelId
									for each (pair<pathNode*, int> pair2 in pixelLevelConnect)
									{
										if (nextId == pair2.second){
											node->directedConnects[pair1.first] = pair2.first;
										}
									
									}

									
								}
							}

						}
					}
				}
				cout << endl;
				default_random_engine e;
				uniform_real_distribution<float> u(0.0f, 1.0f);
				
				for each (pair<pathNode*, pathNode*> direct in startDirects)
				{
					int pathLength = 0;
					pathNode* sNode = direct.first;
					for (int n = 0; n < sNode->connects.size(); n++)
					{
						pathNode *nowNode = sNode->connects[n];
						if (nowNode != direct.second){
							continue;
						}
						
						Vector4f pathColor(u(e), u(e), u(e), 1.0f);
						pathNode* lastNode = sNode;
						//只畫一次的從lastNode到nowNode
						if (lastNode->path[n].size()>0){
							showEdges.push_back(toVector3f(lastNode->position));
							axisColor.push_back(colorTable[3]);
							sequenceColors.push_back(pathColor);
							showEdges.push_back(toVector3f(lastNode->path[n].front()));
							axisColor.push_back(colorTable[3]);
							sequenceColors.push_back(pathColor);
							for (int i = 0; i < lastNode->path[n].size() - 1; i++){
								showEdges.push_back(toVector3f(lastNode->path[n][i]));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
								showEdges.push_back(toVector3f(lastNode->path[n][i + 1]));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
							}
							showEdges.push_back(toVector3f(lastNode->path[n].back()));
							axisColor.push_back(colorTable[3]);
							sequenceColors.push_back(pathColor);
							showEdges.push_back(toVector3f(nowNode->position));
							axisColor.push_back(colorTable[3]);
							sequenceColors.push_back(pathColor);
						}
						else{
							showEdges.push_back(toVector3f(lastNode->position));
							axisColor.push_back(colorTable[3]);
							sequenceColors.push_back(pathColor);
							showEdges.push_back(toVector3f(nowNode->position));
							axisColor.push_back(colorTable[3]);
							sequenceColors.push_back(pathColor);
						}
						
						pathLength++;

						while (nowNode->directedConnects[lastNode] != NULL)
						{
							
							pathNode* nextNode = nowNode->directedConnects[lastNode];
							int nextIndex;
							for (nextIndex = 0; nextIndex < nowNode->connects.size(); nextIndex++){
								if (nowNode->connects[nextIndex] == nextNode){
									break;
								}
							}
							if (nowNode->path[nextIndex].size()>0){
								showEdges.push_back(toVector3f(nowNode->position));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
								showEdges.push_back(toVector3f(nowNode->path[nextIndex].front()));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
								for (size_t i = 0; i < nowNode->path[nextIndex].size() - 1; i++)
								{
									showEdges.push_back(toVector3f(nowNode->path[nextIndex][i]));
									axisColor.push_back(colorTable[3]);
									sequenceColors.push_back(pathColor);
									showEdges.push_back(toVector3f(nowNode->path[nextIndex][i + 1]));
									axisColor.push_back(colorTable[3]);
									sequenceColors.push_back(pathColor);
								}
								showEdges.push_back(toVector3f(nowNode->path[nextIndex].back()));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
								showEdges.push_back(toVector3f(nextNode->position));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
							}
							else
							{
								showEdges.push_back(toVector3f(nowNode->position));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
								showEdges.push_back(toVector3f(nextNode->position));
								axisColor.push_back(colorTable[3]);
								sequenceColors.push_back(pathColor);
							}
							lastNode = nowNode;
							nowNode = nextNode;
							pathLength++;
						}
						if (nowNode->y_2d == 9){
							cout << "結束于(" << nowNode->x_2d << "," << nowNode->y_2d << ")";
						}
					}
				}
				
			}
		}
		for each(contour* c in *pMap->contours){
			if (c->nextContours.size()>1){
				for (int i = 0; i < c->pfts.size(); i++){
					if (c->id == 30){
						contourColor.push_back(colorTable[2]);
						contourColor.push_back(colorTable[2]);
					}
					else{
						contourColor.push_back(colorTable[0]);
						contourColor.push_back(colorTable[0]);
					}
					float a1 = c->pfts[i].angleByCal;
					contourColor_angle.push_back(Vector4f(a1, a1, a1, 1));
					float a2 = c->pfts[(i + 1) % c->pfts.size()].angleByCal;
					contourColor_angle.push_back(Vector4f(a2, a2, a2, 1));
				}
				for (int i = 0; i < c->pfts.size(); i++){
					divergenceConnects.push_back(toVector3f(c->pfts[i].pointPos));
					divergenceConnects.push_back(toVector3f(c->closetNextPoint[i]));
					int index = 0;
					//cout << "i=" << i << " pointsBelong:" << c->pointsBelong[i];
					for (int idx = 0; idx < c->nextContours.size(); idx++){
						if (c->pointsBelong[i] == c->nextContours[idx]->id)
						{
							connectColors.push_back(colorTable[idx]);
							connectColors.push_back(colorTable[idx]);
							break;
						}
					}
				}
				//cout << "divergenceConnects.size:" << divergenceConnects.size()<<endl;
			}
			else{
				for (int i = 0; i < c->pfts.size(); i++){
					contourColor.push_back(colorTable[0]);
					contourColor.push_back(colorTable[0]);
					float a1 = c->pfts[i].angleByCal;
					contourColor_angle.push_back(Vector4f(a1, a1, a1, 1));
					float a2 = c->pfts[(i + 1) % c->pfts.size()].angleByCal;
					contourColor_angle.push_back(Vector4f(a2, a2, a2, 1));
				}
			}
		}
		cout << "畫出nGraph結束";
		if (createMesh)
		{
			meshCreater mCreater;
			mCreater.create(nGraph, "output.obj");
			cout << "寫入mesh結束";
		}
		method_end = clock();
		cout << "方法耗時:" << ((double)(method_end - method_start)) / CLOCKS_PER_SEC<<endl;
		vector<Vector3f> totalAxisBuffer;
		totalAxisBuffer.insert(totalAxisBuffer.end(),showEdges.begin(),showEdges.end());
		totalAxisBuffer.insert(totalAxisBuffer.end(), showGrids.begin(), showGrids.end());
		totalAxisBuffer.insert(totalAxisBuffer.end(), zeroAxis.begin(), zeroAxis.end());
		totalAxisBuffer.insert(totalAxisBuffer.end(), divergenceConnects.begin(),divergenceConnects.end());
		vector<Vector4f> totalAxisColorBuffer;
		totalAxisColorBuffer.insert(totalAxisColorBuffer.end(), axisColor.begin(), axisColor.end());
		totalAxisColorBuffer.insert(totalAxisColorBuffer.end(), gridColors.begin(), gridColors.end());
		totalAxisColorBuffer.insert(totalAxisColorBuffer.end(), zeroAxisColors.begin(), zeroAxisColors.end());
		totalAxisColorBuffer.insert(totalAxisColorBuffer.end(),connectColors.begin(),connectColors.end());
		//pass1-------------------------------------------------------------------------------------
		glGenVertexArrays(1, VAO);
		cout << "VAO:" << VAO[0] << endl;
		glBindVertexArray(VAO[0]);

		glGenBuffers(2, vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER,Vertices.size()*sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(float), &distPercentage[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);
		shaderHandle = CompileShaders("..\\shader\\shader3.vs", "..\\shader\\shader3.fs");
		
		//pass2------------------------------------------------------------------------------------
		glGenVertexArrays(1, &e_VAO);
		glBindVertexArray(e_VAO);

		glGenBuffers(1, &e_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 6 * sizeof(Vector3f), &edge_Vertex[0], GL_STATIC_DRAW);
		//glBufferData(GL_ARRAY_BUFFER, showEdges.size()*sizeof(Vector3f),&showEdges[0], GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		e_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_black.fs");
		//pass3------------------------------------------------------------------------------------
		glGenVertexArrays(1, &v_VAO);
		glBindVertexArray(v_VAO);

		glGenBuffers(1, &v_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
		glBufferData(GL_ARRAY_BUFFER, showPoints.size()*sizeof(Vector3f), showPoints.data(), GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		v_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_red.fs");
		p_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_blue.fs");
		//pass4------------------------------------------------------------------------------------
		glGenVertexArrays(1, &axis_VAO);
		glBindVertexArray(axis_VAO);

		glGenBuffers(2, axis_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, axis_vbo[0]);
		glBufferData(GL_ARRAY_BUFFER, totalAxisBuffer.size()*sizeof(Vector3f), &totalAxisBuffer[0], GL_STATIC_DRAW); //showEdges.size()*sizeof(Vector3f), &showEdges[0], GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, axis_vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		
		glBindBuffer(GL_ARRAY_BUFFER, axis_vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, totalAxisColorBuffer.size()*sizeof(Vector4f), &totalAxisColorBuffer[0], GL_STATIC_DRAW);//axisColor.size()*sizeof(Vector4f), &axisColor[0], GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, axis_vbo[1]);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
		axis_shaderHandle = CompileShaders("..\\shader\\just_gWorld_Color.vs", "..\\shader\\just_color.fs"); //"..\\shader\\just_green.fs");


		vector<int> shape;
		Vector3f* Varray= toVector3Arrays(lines, shape);
		c_totalNum = 0;
		for (int i = 0; i < shape.size(); i++){
			c_totalNum += shape[i];
		}
		//pass5
		glGenVertexArrays(1, &c_VAO);
		glBindVertexArray(c_VAO);

		glGenBuffers(2, c_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[0]);
		glBufferData(GL_ARRAY_BUFFER, c_totalNum*sizeof(Vector3f), Varray, GL_STATIC_DRAW);
		//cout << "c_totalNum:" << c_totalNum << "contourColor size:" << contourColor.size();
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, contourColor.size()*sizeof(Vector4f), contourColor.data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[1]);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);

	}
	void showGUI(){
		ImGui_ImplOpenGL2_NewFrame();
		ImGui_ImplGLUT_NewFrame();
		if (show_demo_window)
			ImGui::ShowDemoWindow(&show_demo_window);

		ImGui::Begin("test imgui");
		ImGui::Text("This is some useful text.");
		ImGui::End();
		ImGui::Render();
		ImGuiIO& io = ImGui::GetIO();
		glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
	}
	void onRender(){
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//showGUI();
		//glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f), rotateAngle, glm::vec3(0, 1, 0));
		glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f), control_rotate_y, glm::vec3(0, 1, 0));
		glm::mat4 rotat_x = glm::rotate(glm::mat4(1.0f), control_rotate_x, glm::vec3(1, 0, 0));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 0, 0, control_cameraZ }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		//rotateAngle += rotateSpeed;

		glBindVertexArray(VAO[0]);
		glUseProgram(shaderHandle);

		GLuint loc = glGetUniformLocation(shaderHandle, "gWorld");
		assert(loc != 0xFFFFFFFF);
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(dynFakeOrthoMat*glm::mat4(1.0f)));
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		if(showFace)
			glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces() * 9);

		//pass2
		glBindVertexArray(e_VAO);
		glUseProgram(e_shaderHandle);
		GLuint loc2 = glGetUniformLocation(e_shaderHandle, "gWorld");
		assert(loc2 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc2, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glLineWidth(2);
		if (showEdge)
			glDrawArrays(GL_LINES, 0, mesh.n_faces() * 18);
		//pass3
		glBindVertexArray(v_VAO);
		glUseProgram(p_shaderHandle);
		GLuint loc3 = glGetUniformLocation(p_shaderHandle, "gWorld");
		assert(loc3 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc3, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glPointSize(7);
		//cout << "points size:" << showPoints.size()<<endl;
		if (showPoint)
			glDrawArrays(GL_POINTS, 0, showPoints.size());
		//pass4
		glBindVertexArray(c_VAO);
		glUseProgram(axis_shaderHandle);//_shaderHandle);
		GLuint loc4 = glGetUniformLocation(axis_shaderHandle, "gWorld");
		assert(loc4 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc4, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glLineWidth(5);
		if (showContour)
			glDrawArrays(GL_LINES, 0, c_totalNum);
		//pass5
		glBindVertexArray(axis_VAO);
		glUseProgram(axis_shaderHandle);
		GLuint loc5 = glGetUniformLocation(axis_shaderHandle, "gWorld");
		assert(loc5 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc5, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glLineWidth(3);
		if (showAxis){
			if (gridOnly){
				glDrawArrays(GL_LINES, showEdges.size(), showGrids.size());
			}
			else if (drawZeroAxis){
				//cout << "繪製開始:" << showEdges.size() + showGrids.size() << "結束于:" << showEdges.size() + showGrids.size() + zeroAxis.size()<<"zeroAxis.size:"<<zeroAxis.size();
				glDrawArrays(GL_LINES, showEdges.size() + showGrids.size(), zeroAxis.size());
				//glDrawArrays(GL_LINES, showEdges.size() + showGrids.size(), showEdges.size()+showGrids.size()+zeroAxis.size());
			}
			else if (drawDivergenceConnect){
				glDrawArrays(GL_LINES, showEdges.size() + showGrids.size() + zeroAxis.size(), divergenceConnects.size());
			}
			else if (demoSequences)
			{
				if (nowStep>showEdges.size())
				{
					glDrawArrays(GL_LINES, 0, showEdges.size());
				}
				else
				{
					//cout << "第" << nowStep << "幀";
					glDrawArrays(GL_LINES, 0, nowStep++);
					//cout << "=>" << nowStep << "幀. ";
				}
			}
			else
			{
				glDrawArrays(GL_LINES, 0, showEdges.size());
			}
		}
	}
	void onKeyDown(unsigned char key, int mx, int my){
		if (key == 'w' || key == 'W'){
			//cout << "鍵w被按下";
			control_rotate_x -= rotateSpeed;
		}
		else if (key == 'a' || key == 'A'){
			//cout << "鍵a被按下";
			control_rotate_y -= rotateSpeed;
		}
		else if (key == 's' || key == 'S'){
			//cout << "鍵s被按下";
			control_rotate_x += rotateSpeed;
		}
		else if (key == 'd' || key == 'D'){
			//cout << "鍵d被按下";
			control_rotate_y += rotateSpeed;
		}
		else if (key == 'q' || key == 'Q'){
			control_cameraZ += cameraSpeed;
		}
		else if (key == 'e' || key == 'E'){
			control_cameraZ -= cameraSpeed;
		}

		if (key == '1'){
			showFace = true;
			showEdge = true;
			showContour = true;
			showAxis = false;
		}
		else if (key == '2'){
			showFace = false;
			showEdge = false;
			showContour = true;
			showAxis = true;
			gridOnly = false;
			drawDivergenceConnect = false;
			drawZeroAxis = true;
			glBindBuffer(GL_ARRAY_BUFFER, c_vbo[1]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, contourColor_angle.size()*sizeof(Vector4f), contourColor_angle.data());
		}
		else if (key=='3')
		{
			showFace = false;
			showEdge = false;
			showContour = true;
			showAxis = true;
			gridOnly = false;
			drawZeroAxis = false;
			drawDivergenceConnect = true;
			demoSequences = false;
		}
		else if (key == '4')
		{
			showFace = true;
			showEdge = false;
			showContour = true;
			showAxis = true;
			gridOnly = true;
			drawZeroAxis = false;
			drawDivergenceConnect = false;
			demoSequences = false;
			glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, mesh.n_faces() * 3 * sizeof(float), whitePercentage);
		}
		else if (key == '5'){
			showFace = true;
			showEdge = false;
			showContour = true;
			showAxis = true;
			gridOnly = false;
			drawZeroAxis = false;
			drawDivergenceConnect = false;
			demoSequences = false;
			glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
			glBufferSubData(GL_ARRAY_BUFFER,0, mesh.n_faces() * 3 * sizeof(float), whitePercentage);
		}
		else if (key == '6'){
			showFace = true;
			showEdge = false;
			showContour = true;
			showAxis = true;
			gridOnly = false;
			drawZeroAxis = false;
			drawDivergenceConnect = false;
			demoSequences = true;
			glBindBuffer(GL_ARRAY_BUFFER, axis_vbo[1]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sequenceColors.size()*sizeof(Vector4f), &sequenceColors[0]);
		}
	}
};
class showTensor:public script
{
private:
	GLuint VAO;
	GLuint e_VAO;
	GLuint v_VAO;
	//polarMap* pMap;
	tensorMap *tMap;
	int center = 0;
public:
	MyMesh  mesh;

	GLuint vbo[2];
	GLuint e_vbo;
	GLuint v_vbo;

	GLuint shaderHandle;
	GLuint e_shaderHandle;
	GLuint v_shaderHandle;
	float nearZ = 1.0f;
	float farZ = 60.0f;
	float z = 0.5;
	float rotateAngle = 0;
	float rotateSpeed = 0.01f;
	vector<Vector3f> showPoints;
	void onInit(){

		cout << "showTensor onInit" << endl;
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		mesh.update_normals();


		if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
		{
			cerr << "Error: Cannot read mesh:" << filename;
			return;
		}
		else
		{
			cout << "load mesh成功" << endl;
		}



		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}



		mesh.request_vertex_texcoords2D();

		//cout << "点一,归一化之前为:" << mesh.point(MyMesh::VertexHandle(0));
		mesh = normalizeMesh(mesh);

		//創建tMap
		tMap = new tensorMap(mesh, center, kernel_node);
		tMap->createMap();
		tMap->createTensor();
		//pMap = new polarMap(mesh, center);
		//pMap->createMap();
		//pMap->showPoints(1000);
		//pMap->debugDist(200);

		polarPoint farthest = tMap->getPointFrom(center);
		for (int i = 1; i < tMap->pointNum(); i++){
			if (farthest.distance < tMap->getPointFrom(i).distance){//如果第i個點距離比farthest遠則更新farthest
				farthest = tMap->getPointFrom(i);
			}
		}
		cout << "farthest distance:" << farthest.distance << endl;
		Vector3f *Vertices = new Vector3f[mesh.n_faces() * 3];

		//Vector3f *Normals = new Vector3f[mesh.n_faces() * 3];
		Vector3f *edge_Vertex = new Vector3f[mesh.n_faces() * 6];
		float *distPercentage = new float[mesh.n_faces() * 3];

		
		for (int f = 0; f < mesh.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			if (f < 10){
				cout << "face:" << f << ":";
			}
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				Vertices[f * 3 + vindex] = Vector3f(pos[0], pos[1], pos[2]);
				
				distPercentage[f * 3 + vindex] = tMap->getPointFrom((*fv_it).idx()).angle / (2 * M_PI);
				if (f < 10){
					//cout << "頂點" << (*fv_it).idx() << ":" << pMap->getPointFrom((*fv_it).idx()).distance << " ";
					cout << "角度" << (*fv_it).idx() << ":" << tMap->getPointFrom((*fv_it).idx()).angle / (2 * M_PI);
				}
				vindex++;
			}
			if (f < 10){
				cout << endl;
			}
			//第二個pass的邊
			int eindex=0;
			for (MyMesh::FaceHalfedgeCWIter fhe_cit = mesh.fh_cwbegin(fh); fhe_cit != mesh.fh_cwend(fh); fhe_cit++){
				float* from_pos = mesh.point(mesh.from_vertex_handle(*fhe_cit)).data();
				float* to_pos = mesh.point(mesh.to_vertex_handle(*fhe_cit)).data();
				//cout << "index:" << f * 6 + eindex * 2 << "from pos:(" << from_pos[0] << "," << from_pos[1] << "," << from_pos[2] << ")";
				edge_Vertex[f * 6 + eindex * 2] = Vector3f(from_pos[0], from_pos[1], from_pos[2]);
				//cout << "index" << f * 6 + eindex * 2 +1<< "to pos:(" << to_pos[0] << "," << to_pos[1] << "," << to_pos[2] << ")" << endl;
				edge_Vertex[f * 6 + eindex * 2 +1] = Vector3f(to_pos[0], to_pos[1], to_pos[2]);
				eindex++;
			}
			//cout << endl;
		}
		showPoints.clear();
		MyMesh::VertexHandle startPoint = MyMesh::VertexHandle( 1423 );
		OpenMesh::Vec3f point = mesh.point(startPoint);
		//showPoints.push_back(Vector3f(point[0], point[1], point[2]));
		for (int idx = 0; idx < mesh.n_vertices(); idx++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(idx);
			OpenMesh::Vec3f point = mesh.point(vh);
			showPoints.push_back(Vector3f(point[0], point[1], point[2]));
		}
		//完成處理mesh

		//完成創建map
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(float), &distPercentage[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);
		shaderHandle = CompileShaders("..\\shader\\shader3.vs", "..\\shader\\shader3.fs");
		//pass2------------------------------------------------------------------------------------
		glGenVertexArrays(1, &e_VAO);
		glBindVertexArray(e_VAO);

		glGenBuffers(1,&e_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 6 * sizeof(Vector3f), &edge_Vertex[0], GL_STATIC_DRAW);
		
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		e_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_black.fs");
		//pass3------------------------------------------------------------------------------------
		glGenVertexArrays(1, &v_VAO);
		glBindVertexArray(v_VAO);

		glGenBuffers(1, &v_vbo);
		glBindBuffer(GL_ARRAY_BUFFER,v_vbo);
		glBufferData(GL_ARRAY_BUFFER, showPoints.size()*sizeof(Vector3f),showPoints.data(),GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		v_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_red.fs");
	}
	void onRender(){
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f), rotateAngle, glm::vec3(1, 0, 0));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 1, 0, z }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		//rotateAngle += rotateSpeed;

		glBindVertexArray(VAO);
		glUseProgram(shaderHandle);
		
		GLuint loc = glGetUniformLocation(shaderHandle, "gWorld");
		assert(loc != 0xFFFFFFFF);
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(dynFakeOrthoMat*glm::mat4(1.0f)));
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y));
		glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces() * 9);
		
	//pass2
		glBindVertexArray(e_VAO);
		glUseProgram(e_shaderHandle);
		GLuint loc2 = glGetUniformLocation(e_shaderHandle, "gWorld");
		assert(loc2 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc2, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y));
		glLineWidth(1);
		glDrawArrays(GL_LINES, 0, mesh.n_faces() * 18);
	//pass3
		glBindVertexArray(v_VAO);
		glUseProgram(v_shaderHandle);
		GLuint loc3 = glGetUniformLocation(v_shaderHandle, "gWorld");
		assert(loc3 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc3, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y));
		glPointSize(3);
		//cout << "points size:" << showPoints.size()<<endl;
		glDrawArrays(GL_POINTS, 0, showPoints.size());
	}
};
class showTextureMapping :public script
{
private:
	GLuint VAO;
	GLuint e_VAO;
	GLuint v_VAO;
	//polarMap* pMap;
	//tensorMap *tMap;
	int center = 0;
	float elemPerLevel = 20.0f;
	float levelNumber = 20.0f;
	vector<OpenMesh::Vec2f> vertex_uv;
public:
	MyMesh  mesh;

	GLuint vbo[2];
	GLuint e_vbo;
	GLuint v_vbo;
	GLuint shaderHandle;
	GLuint e_shaderHandle;
	GLuint v_shaderHandle;
	float nearZ = 1.0f;
	float farZ = 60.0f;
	float z = 2;//0.5;
	float rotateAngle = 0;
	float rotateSpeed = 0.005f;
	vector<Vector3f> showPoints;
	GLuint textureID;
	void onInit(){

		cout << "showTensor onInit" << endl;
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();



		if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
		{
			cerr << "Error: Cannot read mesh:" << filename;
			return;
		}
		else
		{
			cout << "load mesh成功" << endl;
		}



		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}
		mesh.update_normals();


		mesh.request_vertex_texcoords2D();

		//cout << "点一,归一化之前为:" << mesh.point(MyMesh::VertexHandle(0));
		mesh = normalizeMesh(mesh);

		OpenMesh::Vec3f meshCenter(0, 0, 0);
		OpenMesh::Vec3f upVertex(0, -9999, 0);
		OpenMesh::Vec3f lowVertex(0, 9999, 0);
		//遍历所有模型顶点加总获得平均值
		for (int i = 0; i < mesh.n_vertices(); i++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
			MyMesh::Point point = mesh.point(vh);
			meshCenter += point;
			if (point[1] > upVertex[1]){
				upVertex = point;
			}
			if (point[1] < lowVertex[1]){
				lowVertex = point;
			}

		}
		meshCenter /= mesh.n_vertices();

		for (int i = 0; i < mesh.n_vertices(); i++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
			OpenMesh::Vec3f point = mesh.point(vh);
			OpenMesh::Vec3f vec = point - meshCenter;
			OpenMesh::Vec2f vec_in_xz = OpenMesh::Vec2f(vec[0], vec[2]).normalize();
			float angle = 0;
			//cout << "i=" << i << ":";
			angle = acosf(vec_in_xz[0]) / (2 * M_PI);
			if (vec_in_xz[1] < 0){//說明角度應該在180-360之間
				if (angle < 0.5){
					angle = 1 - angle;
				}
				//cout << "0-180:";
			}
			else{//角度在0-180之間
				//啥也不用做.
			}
			//cout << "angle:" << (angle / (2 * M_PI)) * 360 << "; ";
			float height = (point[1] - lowVertex[1]) / (upVertex[1] - lowVertex[1]);
			vertex_uv.push_back(OpenMesh::Vec2f(angle*elemPerLevel,height*levelNumber));
		}

		Vector3f *Vertices = new Vector3f[mesh.n_faces() * 3];
		Vector2f *uvCoord = new Vector2f[mesh.n_faces() * 3];

		//Vector3f *Normals = new Vector3f[mesh.n_faces() * 3];
		Vector3f *edge_Vertex = new Vector3f[mesh.n_faces() * 6];
		//float *distPercentage = new float[mesh.n_faces() * 3];


		for (int f = 0; f < mesh.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			if (f < 10){
				cout << "face:" << f << ":";
			}
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				int vidx = (*fv_it).idx();
				Vertices[f * 3 + vindex] = Vector3f(pos[0], pos[1], pos[2]);
				uvCoord[f * 3 + vindex] = Vector2f(vertex_uv[vidx][0], vertex_uv[vidx][1]);
				vindex++;
			}
			if (f < 10){
				cout << endl;
			}
			//第二個pass的邊
			int eindex = 0;
			for (MyMesh::FaceHalfedgeCWIter fhe_cit = mesh.fh_cwbegin(fh); fhe_cit != mesh.fh_cwend(fh); fhe_cit++){
				float* from_pos = mesh.point(mesh.from_vertex_handle(*fhe_cit)).data();
				float* to_pos = mesh.point(mesh.to_vertex_handle(*fhe_cit)).data();
				//cout << "index:" << f * 6 + eindex * 2 << "from pos:(" << from_pos[0] << "," << from_pos[1] << "," << from_pos[2] << ")";
				edge_Vertex[f * 6 + eindex * 2] = Vector3f(from_pos[0], from_pos[1], from_pos[2]);
				//cout << "index" << f * 6 + eindex * 2 +1<< "to pos:(" << to_pos[0] << "," << to_pos[1] << "," << to_pos[2] << ")" << endl;
				edge_Vertex[f * 6 + eindex * 2 + 1] = Vector3f(to_pos[0], to_pos[1], to_pos[2]);
				eindex++;
			}
			//cout << endl;
		}
		showPoints.clear();
		MyMesh::VertexHandle startPoint = MyMesh::VertexHandle(1423);
		OpenMesh::Vec3f point = mesh.point(startPoint);
		for (int idx = 0; idx < mesh.n_vertices(); idx++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(idx);
			OpenMesh::Vec3f point = mesh.point(vh);
			showPoints.push_back(Vector3f(point[0], point[1], point[2]));
		}
		//貼圖
		string texturePath = "../Resources/texture_for_elem1.png";//"../Resources/voidTexture.jpg";//"../Resources/elem1.png";
		cv::Mat image = cv::imread(texturePath);
		cout << "Image size:(" << image.size().width << "," << image.size().height << ")";
		cv::flip(image, image, 0);
		cv::imshow("Texture", image);
		cv::waitKey(0);
		cv::destroyAllWindows();
		glGenTextures(1,&textureID);
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image.size().width, image.size().height, 0, GL_BGR, GL_UNSIGNED_BYTE, image.data);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);



		//完成創建map
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector2f), &uvCoord[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
		shaderHandle = CompileShaders("..\\shader\\withTexture.vs", "..\\shader\\withTexture.fs");
		//pass2------------------------------------------------------------------------------------
		glGenVertexArrays(1, &e_VAO);
		glBindVertexArray(e_VAO);

		glGenBuffers(1, &e_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 6 * sizeof(Vector3f), &edge_Vertex[0], GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		e_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_black.fs");
		//pass3------------------------------------------------------------------------------------
		glGenVertexArrays(1, &v_VAO);
		glBindVertexArray(v_VAO);

		glGenBuffers(1, &v_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
		glBufferData(GL_ARRAY_BUFFER, showPoints.size()*sizeof(Vector3f), showPoints.data(), GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		v_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_red.fs");
	}
	void onRender(){
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f), rotateAngle, glm::vec3(0, 1, 0));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 1, 0, z }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });
		rotateAngle += rotateSpeed;

		glBindVertexArray(VAO);
		glUseProgram(shaderHandle);

		GLuint loc = glGetUniformLocation(shaderHandle, "gWorld");
		assert(loc != 0xFFFFFFFF);
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(dynFakeOrthoMat*glm::mat4(1.0f)));
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y));
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, textureID);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		GLuint textureLoc = glGetUniformLocation(shaderHandle, "myTextureSampler");
		glUniform1i(textureLoc, 0);
		glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces() * 9);

		//pass2
		glBindVertexArray(e_VAO);
		glUseProgram(e_shaderHandle);
		GLuint loc2 = glGetUniformLocation(e_shaderHandle, "gWorld");
		assert(loc2 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc2, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y));
		glLineWidth(1);
		glDrawArrays(GL_LINES, 0, mesh.n_faces() * 18);
		//pass3
		glBindVertexArray(v_VAO);
		glUseProgram(v_shaderHandle);
		GLuint loc3 = glGetUniformLocation(v_shaderHandle, "gWorld");
		assert(loc3 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc3, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y));
		glPointSize(3);
		//cout << "points size:" << showPoints.size()<<endl;
		//glDrawArrays(GL_POINTS, 0, showPoints.size());
	}
};
class showProjectionPath :public script
{
private:
	GLuint VAO;
	GLuint e_VAO;
	GLuint v_VAO;
	//polarMap* pMap;
	//tensorMap *tMap;
	int center = 0;
	float elemPerLevel = 20.0f;
	float levelNumber = 20.0f;
	vector<OpenMesh::Vec2f> vertex_uv;
public:
	MyMesh  mesh;

	GLuint vbo[2];
	GLuint e_vbo;
	GLuint v_vbo[2];
	GLuint shaderHandle;
	GLuint e_shaderHandle;
	GLuint v_shaderHandle;
	float nearZ = 1.0f;
	float farZ = 60.0f;
	float z = 2;//0.5;

	vector<Vector3f> showPoints;
	vector<Vector4f> pathColors;
	GLuint textureID;
	coverAxis lastCover;
	const float getPathYAsZeroThreshold = 0.00001f;
	const float matchEndPointThreshold = 0.001f;
	float pathLengthThreshold = 10.0f;
	const int startFace_idx = 247;
	const int endFace_idx = 548;
	const char* modelName = "../Resources/bunny.obj";//"../Resources/horse230630.obj";
protected:
	float control_rotate_y = 0;
	float control_rotate_x = 0;
	float rotateSpeed = 0.01f;//0.0002f;
	float control_cameraZ = 2.0f;
	float cameraSpeed = 0.01f;
	glm::vec3 toGlmVec3(OpenMesh::Vec3f v){
		return glm::vec3(v[0], v[1], v[2]);
	}
	glm::vec3 centerOf(MyMesh::FaceHandle face){
		OpenMesh::Vec3f center(0, 0, 0);
		int Vnum = 0;
		for (MyMesh::FaceVertexCWIter fv_cwit = mesh.fv_cwbegin(face); fv_cwit != mesh.fv_cwend(face);fv_cwit++)
		{
			center += mesh.point(*fv_cwit);
			Vnum++;
		}
		center /= Vnum;
		return toGlmVec3(center);
	}
	vector<OpenMesh::Vec3f> getPointsBtw(int faceStartIdx,int faceEndIdx ,bool &fail){
		clock_t start_clock,end_clock;
		start_clock= clock();
		int stepCount = 0;
		//cout << "startClock:" << start_clock << "; ";
		vector<OpenMesh::Vec3f> points;

		
		MyMesh::FaceHandle face_start = MyMesh::FaceHandle(faceStartIdx);
		glm::vec3 startPoint = centerOf(face_start);
		glm::vec3 norm_start = glm::normalize(toGlmVec3(mesh.normal(face_start)));
		
		MyMesh::FaceHandle face_end = MyMesh::FaceHandle(faceEndIdx);
		glm::vec3 endPoint = centerOf(face_end);
		glm::vec3 norm_end = glm::normalize(toGlmVec3(mesh.normal(face_end)));

		int endFaceIdx = faceEndIdx;
		//cout << "終點fidx:" << endFaceIdx;
		points.push_back(OpenMesh::Vec3f(startPoint[0], startPoint[1], startPoint[2]));
		//cout << "終點為:(" << endPoint[0] << "," << endPoint[1] << "," << endPoint[2] << ") 終點面:" << face_end.idx()<<endl;
		//cout << "getPointsBtw " << contour1.id << " angle:" << at_angle << ":-------------------------------------------" << endl;

		glm::vec3 axis_x = glm::normalize(endPoint - startPoint);

		glm::vec3 temp_z = glm::normalize(norm_start + norm_end);
		//cout << "temp_z:(" << temp_z[0] << "," << temp_z[1] << "," << temp_z[2] << ")" << endl;
		glm::vec3 axis_y = glm::normalize(glm::cross(axis_x, temp_z));
		//cout << "axis_y:(" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;
		glm::vec3 axis_z = glm::normalize(glm::cross(axis_y, axis_x));
		//cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;

		coverAxis trans = coverAxis(startPoint, axis_x, axis_y, axis_z);
		MyMesh::FaceHandle nowFace = MyMesh::FaceHandle(faceStartIdx);
		//cout << " 起點fidx:" << nowFace.idx();
		if (nowFace.idx() == endFaceIdx){//如果起點和終點在同一個面上,就不需要找了直接回傳起終點就好了
			points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
			return points;
		}
		glm::vec3 pe = trans.cover(endPoint);
		lastCover = trans;

		MyMesh::HalfedgeHandle mostPossible_he;
		OpenMesh::Vec3f most_interaction;
		float  most_avgX = -9999;
		//float coverDomain_y[2];
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
			MyMesh::VertexHandle from = mesh.from_vertex_handle(*fh_it);
			MyMesh::VertexHandle to = mesh.to_vertex_handle(*fh_it);
			glm::vec3 point1 = trans.cover(toGlmVec3(mesh.point(from)));
			if (abs(point1[1])<getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
			{
				point1[1] = 0;//為了讓演算法能運作,強制設置為0
			}
			glm::vec3 point2 = trans.cover(toGlmVec3(mesh.point(to)));
			if (abs(point2[1])<getPathYAsZeroThreshold)
			{
				point2[1] = 0;//為了讓演算法能運作,強制設置為0
			}
			bool isInteraction = (point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0);//如果有相交的話兩個點會有一個端點y值為正另一個端點y值為負
			if (isInteraction){
				glm::vec3 offset = point2 - point1;
				float time = (0 - point1[1]) / offset[1];
				glm::vec3 coverDomain_intersection = point1 + time*offset;//轉換空間的交點
				float avgX = coverDomain_intersection[0];//轉換空間交點的x值
				if (avgX>most_avgX){//我們取x最遠的轉換空間交點做為第一個交點
					mostPossible_he = *fh_it;
					most_avgX = avgX;
					OpenMesh::Vec3f interaction = mesh.point(from) + (mesh.point(to) - mesh.point(from))*time;
					most_interaction = interaction;
				}
			}
		}
		//cout << "搜尋點開始:" << endl;
		points.push_back(most_interaction);
		MyMesh::HalfedgeHandle now_he = mesh.opposite_halfedge_handle(mostPossible_he);

		nowFace = mesh.opposite_face_handle(mostPossible_he);


		float totalPathLength = 0;
		while (nowFace.idx() != endFaceIdx){
			MyMesh::VertexHandle now_he_p1 = mesh.from_vertex_handle(now_he);

			glm::vec3 now_he_p1_cover = trans.cover(toGlmVec3(mesh.point(now_he_p1)));
			if (abs(now_he_p1_cover[1]) < getPathYAsZeroThreshold)
			{
				now_he_p1_cover[1] = 0;//為了讓演算法能運作,強制設置為0
			}
			MyMesh::VertexHandle now_he_p2 = mesh.to_vertex_handle(now_he);
			glm::vec3 now_he_p2_cover = trans.cover(toGlmVec3(mesh.point(now_he_p2)));
			if (abs(now_he_p2_cover[1]) < getPathYAsZeroThreshold)
			{
				now_he_p2_cover[1] = 0;//為了讓演算法能運作,強制設置為0
			}
			OpenMesh::Vec3f lastPoint = points.back();
			glm::vec3 lastPointCovered = trans.cover(toGlmVec3(lastPoint));
			if (now_he_p1_cover[1] == 0)//如果上一個點是某個halfEdge的端點的話,就不能走普通的路線,因為當交點位於端點時,以這個點為中心向外的任何一條半邊都會滿足條件(因為其中一個點y==0,另一個點不論大於還是小於都會成立),最終導致交點一直卡在那個端點上
			{
				MyMesh::VertexHandle now_vh = now_he_p1;
				vector<MyMesh::VertexHandle> vNeighbors;
				for (MyMesh::VertexVertexCWIter vvcw_it = mesh.vv_cwbegin(now_vh); vvcw_it != mesh.vv_cwend(now_vh); vvcw_it++)
				{
					vNeighbors.push_back(*vvcw_it);
				}
				float max_x = -9999;
				MyMesh::HalfedgeHandle mostPossible_he;
				OpenMesh::Vec3f most_intersection;
				//glm::vec3 most_coverDomain_intersection;
				for (size_t i = 0; i < vNeighbors.size(); i++)
				{
					OpenMesh::Vec3f from = mesh.point(vNeighbors[i]);
					OpenMesh::Vec3f to = mesh.point(vNeighbors[(i + 1) % vNeighbors.size()]);

					glm::vec3 point1 = trans.cover(toGlmVec3(from));
					if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
					{
						point1[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
					glm::vec3 point2 = trans.cover(toGlmVec3(to));
					if (abs(point2[1]) < getPathYAsZeroThreshold)
					{
						point2[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//if (nowFace.idx() == 7297)
					//cout << "vNeighbors" << vNeighbors[i] << ":(" << point1[0] << "," << point1[1] << ") vNeighbors" << vNeighbors[(i + 1) % vNeighbors.size()] << ":(" << point2[0] << "," << point2[1] << ")";
					if ((point1[1] <= 0 && point2[1] >= 0) || (point1[1] >= 0 && point2[1] <= 0))
					{
						glm::vec3 P1toP2 = point2 - point1;
						float time = (0 - point1[1]) / P1toP2[1];
						glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
						//if (nowFace.idx() == 7297)
						//cout << "coverd:(" << coverDomain_intersection[0] << "," << coverDomain_intersection[1] << ")";
						//glm::vec3 interaction = point1 + P1toP2*time;
						if (coverDomain_intersection[0]>max_x)
						{
							//if (nowFace.idx() == 7297)
							//cout << " set!";
							most_intersection = from + (to - from)*time;
							max_x = coverDomain_intersection[0];
							mostPossible_he = mesh.find_halfedge(vNeighbors[i], vNeighbors[(i + 1) % vNeighbors.size()]);
						}
					}
					//cout << endl;
				}

				if (max_x>pe[0])//如果延伸距離超過終點所在的距離
				{
					float time = (pe[0] - lastPointCovered[0]) / (max_x - lastPointCovered[0]);
					if (time < 0){//早已超過endPoint
						fail = true;
						return points;
					}
					OpenMesh::Vec3f possibleEndPoint = time*(most_intersection - lastPoint) + lastPoint;
					if (glm::length(toGlmVec3(possibleEndPoint) - endPoint)<matchEndPointThreshold){//如果可能的交點非常靠近終點了,就認為它是終點
						points.push_back(possibleEndPoint);
						break;
					}
				}
				now_he = mostPossible_he;//這裡不用反向
				nowFace = mesh.face_handle(now_he);
				points.push_back(most_intersection);
				//cout << "over" << endl;
			}
			else if (now_he_p2_cover[1] == 0){
				MyMesh::VertexHandle now_vh = now_he_p2;
				vector<MyMesh::VertexHandle> vNeighbors;
				for (MyMesh::VertexVertexCWIter vvcw_it = mesh.vv_cwbegin(now_vh); vvcw_it != mesh.vv_cwend(now_vh); vvcw_it++)
				{
					vNeighbors.push_back(*vvcw_it);
				}
				float max_x = -9999;
				MyMesh::HalfedgeHandle mostPossible_he;
				OpenMesh::Vec3f most_intersection;
				//glm::vec3 most_coverDomain_intersection;
				for (size_t i = 0; i < vNeighbors.size(); i++)
				{
					OpenMesh::Vec3f from = mesh.point(vNeighbors[i]);
					OpenMesh::Vec3f to = mesh.point(vNeighbors[(i + 1) % vNeighbors.size()]);

					glm::vec3 point1 = trans.cover(toGlmVec3(from));
					if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
					{
						point1[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
					glm::vec3 point2 = trans.cover(toGlmVec3(to));
					if (abs(point2[1]) < getPathYAsZeroThreshold)
					{
						point2[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//if (nowFace.idx() == 7297)
					//cout << "vNeighbors" << vNeighbors[i] << ":(" << point1[0] << "," << point1[1] << ") vNeighbors" << vNeighbors[(i + 1) % vNeighbors.size()] << ":(" << point2[0] << "," << point2[1] << ")";
					if ((point1[1] <= 0 && point2[1] >= 0) || (point1[1] >= 0 && point2[1] <= 0))
					{
						glm::vec3 P1toP2 = point2 - point1;
						float time = (0 - point1[1]) / P1toP2[1];
						glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
						//if (nowFace.idx() == 7297)
						//cout << "coverd:(" << coverDomain_intersection[0] << "," << coverDomain_intersection[1] << ")";
						//glm::vec3 interaction = point1 + P1toP2*time;
						if (coverDomain_intersection[0]>max_x)
						{
							//if (nowFace.idx() == 7297)
							//cout << " set!";
							most_intersection = from + (to - from)*time;
							max_x = coverDomain_intersection[0];
							mostPossible_he = mesh.find_halfedge(vNeighbors[i], vNeighbors[(i + 1) % vNeighbors.size()]);
						}
					}
					//cout << endl;
				}

				if (max_x>pe[0])//如果延伸距離超過終點所在的距離
				{
					float time = (pe[0] - lastPointCovered[0]) / (max_x - lastPointCovered[0]);
					if (time < 0){//早已超過endPoint
						fail = true;
						return points;
					}
					OpenMesh::Vec3f possibleEndPoint = time*(most_intersection - lastPoint) + lastPoint;
					if (glm::length(toGlmVec3(possibleEndPoint) - endPoint)<matchEndPointThreshold){//如果可能的交點非常靠近終點了,就認為它是終點
						points.push_back(possibleEndPoint);
						break;
					}
				}
				now_he = mostPossible_he;//這裡不用反向
				nowFace = mesh.face_handle(now_he);
				points.push_back(most_intersection);
				//cout << "over" << endl;
			}
			else
			{
				for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){

					if ((*fh_it) != now_he){//忽略更新的來源邊
						//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
						OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
						OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
						glm::vec3 point1 = trans.cover(toGlmVec3(from));
						if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
						{
							point1[1] = 0;//為了讓演算法能運作,強制設置為0
						}
						//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
						glm::vec3 point2 = trans.cover(toGlmVec3(to));
						if (abs(point2[1]) < getPathYAsZeroThreshold)
						{
							point2[1] = 0;//為了讓演算法能運作,強制設置為0
						}
						//if (contour1_index == 7)
						//	cout << "interaction v1:" << mesh.from_vertex_handle(*fh_it) << ">(" << point1[0] << "," << point1[1] << ") v2:" << mesh.to_vertex_handle(*fh_it) << ">(" << point2[0] << "," << point2[1] << ")";

						//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";

						if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0)){//找到一個符合兩端的各在x軸一邊的邊
							//cout << "從半邊:" << (*fh_it).idx();
							glm::vec3 P1toP2 = point2 - point1;
							float time = (0 - point1[1]) / P1toP2[1];
							//glm::vec3 interaction = point1 + P1toP2*time;
							OpenMesh::Vec3f intersection = from + (to - from)*time;
							//glm::vec3 intersection_trans = point1 + (point2 - point1)*time;
							totalPathLength += (intersection - points.back()).length();

							points.push_back(intersection);
							//cout << "(" << intersection_trans[0] << "," << intersection_trans[1] << "," << intersection_trans[2] << ")";

							if (totalPathLength > pathLengthThreshold){//說明超過了終點所在的範圍,在雙層結構下很可能會發生,終點在上層,路徑卻走下層
								fail = true;
								return points;
							}
							//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";

							now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
							//if (contour1_index == 7)
							//	cout << "更新半邊到:" << now_he.idx() << endl;
							nowFace = mesh.face_handle(now_he);//推進面
							//coverDomain_y[0] = point1[1];
							//coverDomain_y[1] = point2[1];
							break;
						}
					}
					//cout << ";     ";
					//cout << endl;
				}
			}
			stepCount++;
		}
		//cout << endl;
		points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
		end_clock = clock();
		//cout << "end_clock:" << end_clock;
		double costTime = ((double)(end_clock - start_clock)) / CLOCKS_PER_SEC;
		cout << "getPointsBtw 耗時:" << costTime<<endl;
		cout << "經過步數:" << stepCount << endl;
		//cout << "getPointsBtw end----------------------------------------------"<<endl;
		fail = false;
		return points;
	}
	void onInit(){

		cout << "showTensor onInit" << endl;
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();



		if (!OpenMesh::IO::read_mesh(mesh, modelName, opt))
		{
			cerr << "Error: Cannot read mesh:" << filename;
			return;
		}
		else
		{
			cout << "load mesh成功" << endl;
		}



		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}
		mesh.update_normals();


		mesh.request_vertex_texcoords2D();

		//cout << "点一,归一化之前为:" << mesh.point(MyMesh::VertexHandle(0));
		mesh = normalizeMesh(mesh);
		OpenMesh::Vec3f meshCenter(0, 0, 0);
		OpenMesh::Vec3f upVertex(0, -9999, 0);
		OpenMesh::Vec3f lowVertex(0, 9999, 0);
		//遍历所有模型顶点加总获得平均值
		for (int i = 0; i < mesh.n_vertices(); i++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
			MyMesh::Point point = mesh.point(vh);
			meshCenter += point;
			if (point[1] > upVertex[1]){
				upVertex = point;
			}
			if (point[1] < lowVertex[1]){
				lowVertex = point;
			}

		}
		meshCenter /= mesh.n_vertices();

		for (int i = 0; i < mesh.n_vertices(); i++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
			OpenMesh::Vec3f point = mesh.point(vh);
			OpenMesh::Vec3f vec = point - meshCenter;
			OpenMesh::Vec2f vec_in_xz = OpenMesh::Vec2f(vec[0], vec[2]).normalize();
			float angle = 0;
			//cout << "i=" << i << ":";
			angle = acosf(vec_in_xz[0]) / (2 * M_PI);
			if (vec_in_xz[1] < 0){//說明角度應該在180-360之間
				if (angle < 0.5){
					angle = 1 - angle;
				}
				//cout << "0-180:";
			}
			else{//角度在0-180之間
				//啥也不用做.
			}
			//cout << "angle:" << (angle / (2 * M_PI)) * 360 << "; ";
			float height = (point[1] - lowVertex[1]) / (upVertex[1] - lowVertex[1]);
			vertex_uv.push_back(OpenMesh::Vec2f(angle*elemPerLevel, height*levelNumber));
		}

		Vector3f *Vertices = new Vector3f[mesh.n_faces() * 3];
		Vector2f *uvCoord = new Vector2f[mesh.n_faces() * 3];

		//Vector3f *Normals = new Vector3f[mesh.n_faces() * 3];
		Vector3f *edge_Vertex = new Vector3f[mesh.n_faces() * 6];
		//float *distPercentage = new float[mesh.n_faces() * 3];


		for (int f = 0; f < mesh.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			if (f < 10){
				cout << "face:" << f << ":";
			}
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				int vidx = (*fv_it).idx();
				Vertices[f * 3 + vindex] = Vector3f(pos[0], pos[1], pos[2]);
				uvCoord[f * 3 + vindex] = Vector2f(vertex_uv[vidx][0], vertex_uv[vidx][1]);
				//MyMesh::Normal norm = mesh.normal(*fv_it);
				//Normals[f * 3 + vindex] = Vector3f(norm[0], norm[1], norm[2]);
				//if ((*fv_it).idx() == 0){
				//	cout<<"0 distance:"<< (pMap->getPointFrom((*fv_it).idx())).distance;
				//}

				//distPercentage[f * 3 + vindex] = pMap->getPointFrom((*fv_it).idx()).distance/farthest.distance;
				//distPercentage[f * 3 + vindex] = 0.75f;//tMap->getPointFrom((*fv_it).idx()).angle / (2 * M_PI);
				vindex++;
			}
			if (f < 10){
				cout << endl;
			}
			/*for (MyMesh::FaceEdgeCWIter fe_cit = mesh.fe_cwbegin(fh); fe_cit != mesh.fe_cwend(fh); fe_cit++){
			MyMesh::EdgeHandle eh = *fe_cit;
			mesh.ev
			}*/
			//第二個pass的邊
			int eindex = 0;
			for (MyMesh::FaceHalfedgeCWIter fhe_cit = mesh.fh_cwbegin(fh); fhe_cit != mesh.fh_cwend(fh); fhe_cit++){
				float* from_pos = mesh.point(mesh.from_vertex_handle(*fhe_cit)).data();
				float* to_pos = mesh.point(mesh.to_vertex_handle(*fhe_cit)).data();
				//cout << "index:" << f * 6 + eindex * 2 << "from pos:(" << from_pos[0] << "," << from_pos[1] << "," << from_pos[2] << ")";
				edge_Vertex[f * 6 + eindex * 2] = Vector3f(from_pos[0], from_pos[1], from_pos[2]);
				//cout << "index" << f * 6 + eindex * 2 +1<< "to pos:(" << to_pos[0] << "," << to_pos[1] << "," << to_pos[2] << ")" << endl;
				edge_Vertex[f * 6 + eindex * 2 + 1] = Vector3f(to_pos[0], to_pos[1], to_pos[2]);
				eindex++;
			}
			//cout << endl;
		}
		//創建projection path
		clock_t start_clock, end_clock;
		bool fail = true;
		start_clock = clock();
		vector<OpenMesh::Vec3f> path = getPointsBtw(startFace_idx, endFace_idx,fail);
		end_clock = clock();
		cout << "getPointsBtw fail:" << fail;
		cout << "投影路徑耗時:" << ((double)(end_clock - start_clock)) / CLOCKS_PER_SEC<<endl;
		for (size_t i = 0;i<0;i++) //i < path.size()-1; i++)
		{
			showPoints.push_back(toVector3f(path[i]));
			pathColors.push_back(Vector4f(1.0f, 0.0f, 0, 1.0f));
			showPoints.push_back(toVector3f(path[i + 1]));
			pathColors.push_back(Vector4f(1.0f, 0.0f, 0, 1.0f));
		}
		
		vector<float> vertexs;
		vector<int> faces;

		for (int v = 0; v < mesh.n_vertices(); v++)
		{
			OpenMesh::Vec3f vpos = mesh.point(OpenMesh::VertexHandle(v));
			vertexs.push_back(vpos[0]);
			vertexs.push_back(vpos[1]);
			vertexs.push_back(vpos[2]);
		}
		for (int fIdx = 0; fIdx < mesh.n_faces(); fIdx++){
			MyMesh::FaceHandle fh = MyMesh::FaceHandle(fIdx);
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it !=mesh.fv_end(fh); fv_it++)
			{
				faces.push_back((*fv_it).idx());
			}
		}

		geodesicPathFinder geodesic(vertexs, faces);
		//start_clock = clock();
		vector<float*> gdc_path= geodesic.getPathBtw_face(startFace_idx, endFace_idx);

		//end_clock = clock();
		//cout << "測地線耗時:" << ((double)(end_clock - start_clock)) / CLOCKS_PER_SEC << endl;
		for (int i = 0; i < gdc_path.size()-1; i++){
			showPoints.push_back(Vector3f(gdc_path[i][0], gdc_path[i][1], gdc_path[i][2]));
			pathColors.push_back(Vector4f(1.0f, 0.0f, 0, 1.0f));
			showPoints.push_back(Vector3f(gdc_path[i+1][0], gdc_path[i+1][1], gdc_path[i+1][2]));
			pathColors.push_back(Vector4f(1.0f, 0.0f, 0, 1.0f));
		}
		string texturePath = "../Resources/voidTexture.jpg";//"../Resources/texture_for_elem1.png";//"../Resources/elem1.png";
		cv::Mat image = cv::imread(texturePath);
		cout << "Image size:(" << image.size().width << "," << image.size().height << ")";
		cv::flip(image, image, 0);
		cv::imshow("Texture", image);
		cv::waitKey(0);
		cv::destroyAllWindows();
		glGenTextures(1, &textureID);
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image.size().width, image.size().height, 0, GL_BGR, GL_UNSIGNED_BYTE, image.data);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		//完成創建map
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector2f), &uvCoord[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
		shaderHandle = CompileShaders("..\\shader\\withTexture.vs", "..\\shader\\withTexture.fs");
		//pass2------------------------------------------------------------------------------------
		glGenVertexArrays(1, &e_VAO);
		glBindVertexArray(e_VAO);

		glGenBuffers(1, &e_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 6 * sizeof(Vector3f), &edge_Vertex[0], GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		e_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_black.fs");
		//pass3------------------------------------------------------------------------------------
		glGenVertexArrays(1, &v_VAO);
		glBindVertexArray(v_VAO);

		glGenBuffers(2, v_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[0]);
		glBufferData(GL_ARRAY_BUFFER, showPoints.size()*sizeof(Vector3f), showPoints.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, pathColors.size()*sizeof(Vector4f), pathColors.data(), GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[1]);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
		v_shaderHandle = CompileShaders("..\\shader\\just_gWorld_Color.vs", "..\\shader\\just_color.fs"); //"..\\shader\\just_green.fs");
	}
	void onRender(){
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f), control_rotate_y, glm::vec3(0, 1, 0));
		glm::mat4 rotat_x = glm::rotate(glm::mat4(1.0f), control_rotate_x, glm::vec3(1, 0, 0));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 0, 0, control_cameraZ }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });


		glBindVertexArray(VAO);
		glUseProgram(shaderHandle);

		GLuint loc = glGetUniformLocation(shaderHandle, "gWorld");
		assert(loc != 0xFFFFFFFF);
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(dynFakeOrthoMat*glm::mat4(1.0f)));
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, textureID);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		GLuint textureLoc = glGetUniformLocation(shaderHandle, "myTextureSampler");
		glUniform1i(textureLoc, 0);
		glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces() * 9);

		//pass2
		glBindVertexArray(e_VAO);
		glUseProgram(e_shaderHandle);
		GLuint loc2 = glGetUniformLocation(e_shaderHandle, "gWorld");
		assert(loc2 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc2, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glLineWidth(1);
		glDrawArrays(GL_LINES, 0, mesh.n_faces() * 18);
		//pass3
		glBindVertexArray(v_VAO);
		glUseProgram(v_shaderHandle);
		GLuint loc3 = glGetUniformLocation(v_shaderHandle, "gWorld");
		assert(loc3 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc3, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		//glPointSize(3);
		//cout << "points size:" << showPoints.size()<<endl;
		glLineWidth(4);
		glDrawArrays(GL_LINES, 0, showPoints.size());
	}
	void onKeyDown(unsigned char key, int mx, int my){
		if (key == 'w' || key == 'W'){
			//cout << "鍵w被按下";
			control_rotate_x -= rotateSpeed;
		}
		else if (key == 'a' || key == 'A'){
			//cout << "鍵a被按下";
			control_rotate_y -= rotateSpeed;
		}
		else if (key == 's' || key == 'S'){
			//cout << "鍵s被按下";
			control_rotate_x += rotateSpeed;
		}
		else if (key == 'd' || key == 'D'){
			//cout << "鍵d被按下";
			control_rotate_y += rotateSpeed;
		}
		else if (key == 'q' || key == 'Q'){
			control_cameraZ += cameraSpeed;
		}
		else if (key == 'e' || key == 'E'){
			control_cameraZ -= cameraSpeed;
		}
	}
};
class IS_Match_Compare :public script
{
private:
	GLuint VAO;
	GLuint e_VAO;
	GLuint v_VAO;
	GLuint c_VAO;
	//polarMap* pMap;
	//tensorMap *tMap;
	int center = 0;//296;
	float elemPerLevel = 20.0f;
	float levelNumber = 20.0f;
	vector<OpenMesh::Vec2f> vertex_uv;
public:
	MyMesh  mesh;

	GLuint vbo[2];
	GLuint e_vbo;
	GLuint v_vbo[2];
	GLuint c_vbo[2];
	GLuint shaderHandle;
	GLuint e_shaderHandle;
	GLuint v_shaderHandle;
	float nearZ = 1.0f;
	float farZ = 60.0f;
	float z = 2;//0.5;

	vector<Vector3f> showPoints;
	vector<Vector4f> pathColors;
	vector<Vector3f> contourPoints;
	vector<Vector4f> contourPoints_color;
	GLuint textureID;
	coverAxis lastCover;
	const char* modelName = "../Resources/bunny.obj";

protected:
	float control_rotate_y = 0;
	float control_rotate_x = 0;
	float rotateSpeed = 0.01f;//0.0002f;
	float control_cameraZ = 2.0f;
	float cameraSpeed = 0.01f;
	glm::vec3 toGlmVec3(OpenMesh::Vec3f v){
		return glm::vec3(v[0], v[1], v[2]);
	}
	glm::vec3 centerOf(MyMesh::FaceHandle face){
		OpenMesh::Vec3f center(0, 0, 0);
		int Vnum = 0;
		for (MyMesh::FaceVertexCWIter fv_cwit = mesh.fv_cwbegin(face); fv_cwit != mesh.fv_cwend(face); fv_cwit++)
		{
			center += mesh.point(*fv_cwit);
			Vnum++;
		}
		center /= Vnum;
		return toGlmVec3(center);
	}
	void showShape(){
		const int length1 = 16;
		float shape1[length1][2] = { { 235, 155 }, { 277, 73 }, { 347, 53 }, { 415, 87 }, { 457, 169 },
		{ 419, 243 }, { 361, 303 }, { 304, 358 }, { 218, 406 }, { 147, 360 },
		{ 95, 306 }, { 42, 246 }, { 19, 170 }, { 65, 95 }, { 115, 60 },
		{ 177, 84 } };
		const int length2 = 14;
		float shape2[length2][2] = { { 235, 155 }, { 277, 73 }, { 347, 53 }, { 415, 87 }, { 457, 169 },
		{ 419, 243 }, { 361, 303 }, { 304, 358 }, { 218, 406 }, {140,407},
		{75, 407}, { 75, 269 }, { 75, 123 }, {140,123} };
		vector<OpenMesh::Vec3f> line1(length1);
		for (int p = 0; p < length1; p++)
		{
			line1[p] = OpenMesh::Vec3f(shape1[p][0], shape1[p][1], 0);
		}
		IS_FeatureMatrix feature1(line1);
		//int feature1_width = feature1.matrix.size[0];
		cv::Mat image1(500 , 500, CV_8UC3, cv::Scalar(0, 0, 0));
		for (size_t i = 0; i < length1; i++)
		{
			float* node = shape1[i];
			cv::circle(image1,cv::Point2i(node[0],node[1]),3,cv::Scalar(0,0,255),2);
			float* node_next = shape1[(i + 1) % length1];
			cv::line(image1, cv::Point2i(node[0], node[1]), cv::Point2i(node_next[0], node_next[1]), cv::Scalar(200, 200, 200), 2);
		}

		
		cv::imshow("shapeImage",image1);
		cv::waitKey(0);
		cv::destroyAllWindows();
		cv::imwrite("../IS_test/graph1.jpg",image1);
		cv::Mat demoImg= feature1.showAsImage(20);
		demoImg *= 255;
		cv::imwrite("../IS_test/feature1.jpg", demoImg);
		
		vector<OpenMesh::Vec3f> line2(length2);
		for (int p = 0; p < length2; p++)
		{
			line2[p] = OpenMesh::Vec3f(shape2[p][0], shape2[p][1], 0);
		}
		IS_FeatureMatrix feature2(line2);
		//int feature2_width = feature2.matrix.size[0];
		cv::Mat image2(500, 500, CV_8UC3, cv::Scalar(0, 0, 0));
		for (size_t i = 0; i < length2; i++)
		{
			float* node = shape2[i];
			cv::circle(image2, cv::Point2i(node[0], node[1]), 3, cv::Scalar(0, 0, 255), 2);
			//if (i + 1 < length2){
			float* node_next = shape2[(i + 1)%length2];
			cv::line(image2, cv::Point2i(node[0], node[1]), cv::Point2i(node_next[0], node_next[1]), cv::Scalar(200, 200, 200), 2);
			//}
		}
		cv::imshow("shapeImage", image2);
		cv::waitKey(0);
		cv::destroyAllWindows();
		cv::imwrite("../IS_test/graph2.jpg", image2);
		cv::Mat demoImg2 = feature2.showAsImage(20);
		demoImg2 *= 255;
		cv::imwrite("../IS_test/feature2.jpg", demoImg2);
		vector<Vector2f> new_shape2;
		for (size_t i = 0; i < length2 - 1; i++)
		{
			new_shape2.push_back(Vector2f(shape2[i][0],shape2[i][1]));
			Vector2f interm{(shape2[i][0] + shape2[i + 1][0]) / 2, (shape2[i][1] + shape2[i + 1][1]) / 2};
			new_shape2.push_back(interm);
		}
		new_shape2.push_back(Vector2f(shape2[length2 - 1][0],shape2[length2-1][1]));
		cv::Mat image3(500, 500, CV_8UC3, cv::Scalar(0, 0, 0));
		for (size_t i = 0; i < new_shape2.size(); i++)
		{
			Vector2f node = new_shape2[i];
			cv::circle(image3, cv::Point2i(node.x, node.y), 3, cv::Scalar(0, 0, 255), 2);
			//if (i + 1 < length2){
			Vector2f node_next = new_shape2[(i + 1) % new_shape2.size()];
			cv::line(image3, cv::Point2i(node.x, node.y), cv::Point2i(node_next.x, node_next.y), cv::Scalar(200, 200, 200), 2);
			//}

		}
		vector<OpenMesh::Vec3f> line3(new_shape2.size());
		for (size_t i = 0; i < new_shape2.size(); i++)
		{
			line3[i] = OpenMesh::Vec3f(new_shape2[i].x, new_shape2[i].y, 0);
		}
		IS_FeatureMatrix feature3(line3);

		cv::imshow("shapeImage", image3);
		cv::waitKey(0);
		cv::destroyAllWindows();
		cv::imwrite("../IS_test/graph3.jpg", image3);
		cv::Mat demoImg3 = feature3.showAsImage(20);
		demoImg3 *= 255;
		cv::imwrite("../IS_test/feature3.jpg", demoImg3);
	}
	void onInit(){
		showShape();
		cout << "showTensor onInit" << endl;
		OpenMesh::IO::Options opt;
		mesh.request_face_normals();
		mesh.request_vertex_normals();



		if (!OpenMesh::IO::read_mesh(mesh, modelName, opt))
		{
			cerr << "Error: Cannot read mesh:" << filename;
			return;
		}
		else
		{
			cout << "load mesh成功" << endl;
		}



		if (mesh.has_face_normals())
		{
			cout << "mesh has face normal!" << endl;
		}
		mesh.update_normals();


		mesh.request_vertex_texcoords2D();

		//cout << "点一,归一化之前为:" << mesh.point(MyMesh::VertexHandle(0));
		mesh = normalizeMesh(mesh);
		OpenMesh::Vec3f meshCenter(0, 0, 0);
		OpenMesh::Vec3f upVertex(0, -9999, 0);
		OpenMesh::Vec3f lowVertex(0, 9999, 0);
		//遍历所有模型顶点加总获得平均值
		for (int i = 0; i < mesh.n_vertices(); i++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
			MyMesh::Point point = mesh.point(vh);
			meshCenter += point;
			if (point[1] > upVertex[1]){
				upVertex = point;
			}
			if (point[1] < lowVertex[1]){
				lowVertex = point;
			}

		}
		meshCenter /= mesh.n_vertices();

		for (int i = 0; i < mesh.n_vertices(); i++){
			MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
			OpenMesh::Vec3f point = mesh.point(vh);
			OpenMesh::Vec3f vec = point - meshCenter;
			OpenMesh::Vec2f vec_in_xz = OpenMesh::Vec2f(vec[0], vec[2]).normalize();
			float angle = 0;
			//cout << "i=" << i << ":";
			angle = acosf(vec_in_xz[0]) / (2 * M_PI);
			if (vec_in_xz[1] < 0){//說明角度應該在180-360之間
				if (angle < 0.5){
					angle = 1 - angle;
				}
				//cout << "0-180:";
			}
			else{//角度在0-180之間
				//啥也不用做.
			}
			//cout << "angle:" << (angle / (2 * M_PI)) * 360 << "; ";
			float height = (point[1] - lowVertex[1]) / (upVertex[1] - lowVertex[1]);
			vertex_uv.push_back(OpenMesh::Vec2f(angle*elemPerLevel, height*levelNumber));
		}

		Vector3f *Vertices = new Vector3f[mesh.n_faces() * 3];
		Vector2f *uvCoord = new Vector2f[mesh.n_faces() * 3];

		//Vector3f *Normals = new Vector3f[mesh.n_faces() * 3];
		Vector3f *edge_Vertex = new Vector3f[mesh.n_faces() * 6];
		//float *distPercentage = new float[mesh.n_faces() * 3];


		for (int f = 0; f < mesh.n_faces(); f++){

			MyMesh::FaceHandle fh = MyMesh::FaceHandle(f);
			int vindex = 0;
			if (f < 10){
				cout << "face:" << f << ":";
			}
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); fv_it++)
			{
				float* pos = mesh.point(*fv_it).data();
				int vidx = (*fv_it).idx();
				Vertices[f * 3 + vindex] = Vector3f(pos[0], pos[1], pos[2]);
				uvCoord[f * 3 + vindex] = Vector2f(vertex_uv[vidx][0], vertex_uv[vidx][1]);
				//MyMesh::Normal norm = mesh.normal(*fv_it);
				//Normals[f * 3 + vindex] = Vector3f(norm[0], norm[1], norm[2]);
				//if ((*fv_it).idx() == 0){
				//	cout<<"0 distance:"<< (pMap->getPointFrom((*fv_it).idx())).distance;
				//}

				//distPercentage[f * 3 + vindex] = pMap->getPointFrom((*fv_it).idx()).distance/farthest.distance;
				//distPercentage[f * 3 + vindex] = 0.75f;//tMap->getPointFrom((*fv_it).idx()).angle / (2 * M_PI);
				vindex++;
			}
			if (f < 10){
				cout << endl;
			}
			/*for (MyMesh::FaceEdgeCWIter fe_cit = mesh.fe_cwbegin(fh); fe_cit != mesh.fe_cwend(fh); fe_cit++){
			MyMesh::EdgeHandle eh = *fe_cit;
			mesh.ev
			}*/
			//第二個pass的邊
			int eindex = 0;
			for (MyMesh::FaceHalfedgeCWIter fhe_cit = mesh.fh_cwbegin(fh); fhe_cit != mesh.fh_cwend(fh); fhe_cit++){
				float* from_pos = mesh.point(mesh.from_vertex_handle(*fhe_cit)).data();
				float* to_pos = mesh.point(mesh.to_vertex_handle(*fhe_cit)).data();
				//cout << "index:" << f * 6 + eindex * 2 << "from pos:(" << from_pos[0] << "," << from_pos[1] << "," << from_pos[2] << ")";
				edge_Vertex[f * 6 + eindex * 2] = Vector3f(from_pos[0], from_pos[1], from_pos[2]);
				//cout << "index" << f * 6 + eindex * 2 +1<< "to pos:(" << to_pos[0] << "," << to_pos[1] << "," << to_pos[2] << ")" << endl;
				edge_Vertex[f * 6 + eindex * 2 + 1] = Vector3f(to_pos[0], to_pos[1], to_pos[2]);
				eindex++;
			}
			//cout << endl;
		}
		//創建projection path
		string texturePath = "../Resources/voidTexture.jpg";//"../Resources/texture_for_elem1.png";//"../Resources/elem1.png";
		cv::Mat image = cv::imread(texturePath);
		cout << "Image size:(" << image.size().width << "," << image.size().height << ")";
		cv::flip(image, image, 0);
		//cv::imshow("Texture", image);
		//cv::waitKey(0);
		//cv::destroyAllWindows();
		glGenTextures(1, &textureID);
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image.size().width, image.size().height, 0, GL_BGR, GL_UNSIGNED_BYTE, image.data);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		//創建等高線----------------------------------------------------------
		int contourNum = 20;
		polarMap* pMap = new polarMap_angleProjection(mesh, center);//new polarMap(mesh, center);
		pMap->meshPath = sphere_path;
		pMap->createMap();
		//pMap->drawMeshToImage();

		polarPoint farthest = pMap->getPointFrom(center);
		for (int i = 1; i < pMap->pointNum(); i++){
			if (farthest.distance < pMap->getPointFrom(i).distance){//如果第i個點距離比farthest遠則更新farthest
				farthest = pMap->getPointFrom(i);
			}
		}
		//cout << "farthest distance:" << farthest.distance << endl;
		float *distPercentage = new float[mesh.n_faces() * 3];
		vector<PointCurve> lines = pMap->getContourLine(contourNum);//getContourLine(contourNum, axis);
		vector<vector<int>> axises = pMap->calBaseAxis();
		Vector4f colorTable[5]{Vector4f(1.0, 0.2, 0.2, 1.0), Vector4f(0.2, 1.0, 0.2, 1.0), Vector4f(0.2, 0.2, 1.0, 1.0),Vector4f(1,1,1,0),Vector4f(0.0, 0.0, 0.0, 1.0)};
		for (int i = 0; i < pMap->contours->size(); i++){
			contour *nowC = pMap->contours->operator[](i);
			vector<OpenMesh::Vec3f> vlist = nowC->toVec3fList();
			if (nowC->nextContours.size()>1){
				for (size_t p = 0; p < nowC->pfts.size() - 1; p++)
				{
					OpenMesh::Vec3f p1 = nowC->pfts[p].pointPos;
					int nextId = nowC->pointsBelong[p];
					int idx = -1;
					for (int i = 0; i < nowC->nextContours.size(); i++){
						if (nowC->nextContours[i]->id == nextId)
						{
							idx = i;
						}
					}
					contourPoints.push_back(toVector3f(p1));
					contourPoints_color.push_back(colorTable[idx]);

					OpenMesh::Vec3f p2 = nowC->pfts[p + 1].pointPos;
					nextId = nowC->pointsBelong[p+1];
					idx = -1;
					for (size_t i = 0; i < nowC->nextContours.size(); i++)
					{
						if (nowC->nextContours[i]->id == nextId){
							idx = i;
						}
					}
					contourPoints.push_back(toVector3f(p2));
					contourPoints_color.push_back(colorTable[idx]);
				}
				//for each (contour*  next in nowC->nextContours)
				for (int n = 0; n < nowC->nextContours.size();n++)
				{
					contour* next = nowC->nextContours[n];

						for (size_t p = 0; p < next->pfts.size() - 1; p++)
						{
							OpenMesh::Vec3f p1 = next->pfts[p].pointPos;
							contourPoints.push_back(toVector3f(p1));
							contourPoints_color.push_back(colorTable[n]);
							OpenMesh::Vec3f p2 = next->pfts[p + 1].pointPos;
							contourPoints.push_back(toVector3f(p2));
							contourPoints_color.push_back(colorTable[n]);
						}
					
				}
				map<int, vector<int_pair>> asignGroups;//索引值是對應的下一條等高線的id,值是一個陣列記錄所有劃分給這個等高線的分斷的起點pfts的索引值和終點pfts的索引值
				contour *current = nowC;
				int currentBelongIdx = -1;
				int_pair record;
				int firstIdx = current->pointsBelong[0];
				for (int v = 0; v < current->pfts.size(); v++){
					if (current->pointsBelong[v] != currentBelongIdx){
						record = int_pair(v);
						currentBelongIdx = current->pointsBelong[v];
						//cout << "v=" << v << "創建新record;   ";
					}
					if (v == current->pfts.size() - 1){
						if (currentBelongIdx == firstIdx){//說明首尾是聯通的
							asignGroups[firstIdx][0].values[0] = record.values[0] - current->pfts.size();
							//cout << "v=" << v << "設置初始值;   ";
						}
					}
					else{
						if (current->pointsBelong[v + 1] != currentBelongIdx){
							record.values[1] = v;
							asignGroups[currentBelongIdx].push_back(record);
							//cout << "v=" << v << "記錄于"<< currentBelongIdx << "[" << record.values[0] << "," << record.values[1] << "];   ";
						}
					}
				}
				for each (contour* next in nowC->nextContours)
				{
					vector<int_pair> segms = asignGroups[next->id];
					int indexInCurrent = 0;
					for (int i = 0; i < nowC->nextContours.size(); i++){
						if (nowC->nextContours[i] == next)
						{
							indexInCurrent = i;
						}
					}
					for each (int_pair segm in segms)
					{
						
						IS_FeatureMatrix feature1(vlist, segm[0], segm[1]);
						cout << "subcurve 從" << segm[0] << "到" << segm[1] << "共" <<feature1.matrix.size[0]<<"個node" << endl;
						feature1.showAsImage(5);
						IS_FeatureMatrix feature2(next->toVec3fList());
						cout << "next contour共" << feature2.matrix.size[0] << "個node";
						feature2.showAsImage(5);
						int M = feature1.matrix.size[0];
						int N = feature2.matrix.size[0];
						
						float minDiff =10000;
						int mint = -1;
						if (M < N){
							
							vector<cv::Mat> Integs = getIntegralImages(feature2.matrix, feature1.matrix);
							for (size_t t = 0; t < N; t++)
							{
								float diff = calAvgDiff(0, t, M, Integs);
								if (minDiff>diff)
								{
									minDiff = diff;
									mint = t;
								}
							}
						}
						int_pair segm_(mint, mint + M - 1);
						cout << "對應的子分段為:" << segm_[0] << "," << segm_[1] << ";";
						int count = 0;
						for (int i = segm[0]; i <= segm[1]; i++){
							int index = i;
							if (index < 0){
								index += current->pfts.size();
							}
							contourPoints.push_back(toVector3f(current->pfts[index].pointPos));
							contourPoints_color.push_back(colorTable[indexInCurrent]);
							index = (mint+count)%next->pfts.size();
							contourPoints.push_back(toVector3f(next->pfts[index].pointPos));
							contourPoints_color.push_back(colorTable[indexInCurrent]);
							count++;
						}

					}
				}
				break;
			}
		}
		/*vector<int> shape;
		Vector3f* Varray = toVector3Arrays(lines, shape);
		int c_totalNum = 0;
		for (int i = 0; i < shape.size(); i++){
			c_totalNum += shape[i];
			//cout << "shape[" << i << "]:" << shape[i]<<" ";
		}*/
		//完成創建map
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(2, vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		cout << "size of Vertex:" << sizeof(Vector3f);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector3f), &Vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 3 * sizeof(Vector2f), &uvCoord[0], GL_STATIC_DRAW);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
		shaderHandle = CompileShaders("..\\shader\\withTexture.vs", "..\\shader\\withTexture.fs");
		//pass2------------------------------------------------------------------------------------
		glGenVertexArrays(1, &e_VAO);
		glBindVertexArray(e_VAO);

		glGenBuffers(1, &e_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glBufferData(GL_ARRAY_BUFFER, mesh.n_faces() * 6 * sizeof(Vector3f), &edge_Vertex[0], GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, e_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		e_shaderHandle = CompileShaders("..\\shader\\just_gWorld.vs", "..\\shader\\just_black.fs");
		//pass3------------------------------------------------------------------------------------
		glGenVertexArrays(1, &v_VAO);
		glBindVertexArray(v_VAO);

		glGenBuffers(2, v_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[0]);
		glBufferData(GL_ARRAY_BUFFER, showPoints.size()*sizeof(Vector3f), showPoints.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, pathColors.size()*sizeof(Vector4f), pathColors.data(), GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, v_vbo[1]);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
		v_shaderHandle = CompileShaders("..\\shader\\just_gWorld_Color.vs", "..\\shader\\just_color.fs"); //"..\\shader\\just_green.fs");
		//pass4------------------------------------------------------------------------------------
		glGenVertexArrays(1, &c_VAO);
		glBindVertexArray(c_VAO);

		glGenBuffers(2, c_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[0]);
		glBufferData(GL_ARRAY_BUFFER, contourPoints.size()*sizeof(Vector3f), contourPoints.data(), GL_STATIC_DRAW);
		//cout << "c_totalNum:" << c_totalNum << "contourColor size:" << contourColor.size();
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, contourPoints_color.size()*sizeof(Vector4f), contourPoints_color.data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, c_vbo[1]);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
	}
	void onRender(){
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.7, 0.8, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 rotat_y = glm::rotate(glm::mat4(1.0f), control_rotate_y, glm::vec3(0, 1, 0));
		glm::mat4 rotat_x = glm::rotate(glm::mat4(1.0f), control_rotate_x, glm::vec3(1, 0, 0));
		glm::mat4 viewMat = glm::perspective(45 / 180.0f*3.1415926f, (float)SCR_WIDTH / (float)SCR_HEIGHT, nearZ, farZ)* glm::lookAt(glm::vec3{ 0, 0, control_cameraZ }, glm::vec3{ 0, 0, 0 }, glm::vec3{ 0, 1, 0 });


		glBindVertexArray(VAO);
		glUseProgram(shaderHandle);

		GLuint loc = glGetUniformLocation(shaderHandle, "gWorld");
		assert(loc != 0xFFFFFFFF);
		//glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(dynFakeOrthoMat*glm::mat4(1.0f)));
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, textureID);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		GLuint textureLoc = glGetUniformLocation(shaderHandle, "myTextureSampler");
		glUniform1i(textureLoc, 0);
		//glDrawArrays(GL_TRIANGLES, 0, mesh.n_faces() * 9);

		//pass2
		glBindVertexArray(e_VAO);
		glUseProgram(e_shaderHandle);
		GLuint loc2 = glGetUniformLocation(e_shaderHandle, "gWorld");
		assert(loc2 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc2, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glLineWidth(1);
		glDrawArrays(GL_LINES, 0, mesh.n_faces() * 18);
		//pass3
		glBindVertexArray(v_VAO);
		glUseProgram(v_shaderHandle);
		GLuint loc3 = glGetUniformLocation(v_shaderHandle, "gWorld");
		assert(loc3 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc3, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		//glPointSize(3);
		//cout << "points size:" << showPoints.size()<<endl;
		glLineWidth(4);
		//glDrawArrays(GL_LINES, 0, showPoints.size());
		glBindVertexArray(c_VAO);
		glUseProgram(v_shaderHandle);
		GLuint loc4 = glGetUniformLocation(v_shaderHandle, "gWorld");
		assert(loc4 != 0xFFFFFFFF);
		glUniformMatrix4fv(loc4, 1, GL_FALSE, glm::value_ptr(viewMat*rotat_y*rotat_x));
		glDrawArrays(GL_LINES, 0, contourPoints.size());
	}
	void onKeyDown(unsigned char key, int mx, int my){
		if (key == 'w' || key == 'W'){
			//cout << "鍵w被按下";
			control_rotate_x -= rotateSpeed;
		}
		else if (key == 'a' || key == 'A'){
			//cout << "鍵a被按下";
			control_rotate_y -= rotateSpeed;
		}
		else if (key == 's' || key == 'S'){
			//cout << "鍵s被按下";
			control_rotate_x += rotateSpeed;
		}
		else if (key == 'd' || key == 'D'){
			//cout << "鍵d被按下";
			control_rotate_y += rotateSpeed;
		}
		else if (key == 'q' || key == 'Q'){
			control_cameraZ += cameraSpeed;
		}
		else if (key == 'e' || key == 'E'){
			control_cameraZ -= cameraSpeed;
		}
	}
};
//腳本部分結束下面是主體---------------------------------------------------------------------------
static void RenderSceneCB()
{

	/*glClear(GL_COLOR_BUFFER_BIT);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glDrawArrays(GL_POINTS, 0, 1);

	glEnableVertexAttribArray(0);*/
	//cout << "onRender!!!!!!" << endl;
	currentScript->onRender();
	glutSwapBuffers();
	glutPostRedisplay();
}

static void nothing(){

}
static void InitializeGlutCallbacks()
{
	glutDisplayFunc(nothing);
	glutIdleFunc(RenderSceneCB);
}
static int getdigit(int num) {
	int digit = 0;
	if (num < 0) {
		return -1;
	}
	else {
		while (num>=10)
		{
			//cout << "num 大小為:" << num;
			num /= 10;
			digit += 1;
		} 
	}
	return digit;
}
static void savebytes(string path,unsigned char* data,int length) {
	/*fstream file;
	file.open("data", ios::in);
	if (!file) {
		cout << "檔案無法開啟\n";
	}
	else
	{
		file << data;
		file.close();
		cout << "檔案儲存完成";
	}*/
	ofstream fout(path);
	if (!fout) {
		cout << "無法寫入檔案\n";
		return;
	}
	for (size_t i = 0; i < length; i++)
	{
		fout << data[i];
	}
	//fout << EOF;
	fout.close();
	cout << "寫入完畢";
}
static char* save_img_path = "profile.png";
void mouseCallBack(int button, int state, int x, int y){
	//cout << "mouse call bask 被觸發 buttom"<<button<<" state"<<state<<" "<<x<<","<<y<<endl;
}
void mouseDragCallBack(int x, int y){
	//cout << "mouse drag call back 被觸發 " << x << "," << y<<" "<<endl;
}
void keyCallback(unsigned char key, int mx, int my){
	/*if (key == 'S'||key == 's'){
		cout << "鍵S被按下"<<endl;
		glReadBuffer(GL_FRONT);
		int i = SCR_WIDTH * 3;
		while (i % 4 != 0)  // 补充数据，直到 i 是的倍数
			++i;
		int PixelDataLength = i * SCR_HEIGHT;
		GLubyte* pPixelData = new GLubyte[PixelDataLength];
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
		glReadPixels(0, 0, SCR_WIDTH, SCR_WIDTH, GL_BGR_EXT, GL_UNSIGNED_BYTE, pPixelData);
		string commend = "..\\ImgWriter\\Imghomework";
		int not_zeros = 0;
		for (int p = 0; p < PixelDataLength; p++) {
			if ((int)pPixelData[p]!=0)
			{
				not_zeros++;
			}
		}
		cout<<"本次不為0的pixel數量為:"<<not_zeros<<endl;
		//cout << "11的位數是" << getdigit(11) << "100的位數是" << getdigit(100) << "0的位數是" << getdigit(0);
		int w_dig = getdigit(SCR_WIDTH);
		int h_dig = getdigit(SCR_HEIGHT);

		char* w_chars = new char[10];
		char* h_chars = new char[10];
		sprintf_s(w_chars, 10, "%d", SCR_WIDTH);
		sprintf_s(h_chars, 10, "%d", SCR_HEIGHT);
		char* l_chars = new char[10];
		sprintf_s(l_chars, 10, "%d", PixelDataLength);
		cout << "width:(" << w_chars << ") height:(" << h_chars << ")"<<endl;
		commend += " ";
		commend += w_chars;
		commend += " ";
		commend += h_chars;
		commend += " ";
		commend += save_img_path;
		commend += " ";
		commend += l_chars;
		
		cout << "commend:" + commend<<endl;
		string path = "data";
		cout << "data path:" <<path;
		savebytes(path,pPixelData, PixelDataLength);
		system(commend.data());
		
	}
	else if (key == 'Z' || key == 'z'){
		((showEdge*)currentScript)->updateThreshold(-0.01f);
	}
	else if (key == 'X' || key == 'x')
	{
		((showEdge*)currentScript)->updateThreshold(0.01f);
	}*/
	//cout << "keyCallback被呼叫!";
	currentScript->onKeyDown(key, mx, my);
}

int main(int argc, char** argv)
{

	currentScript =new showPolar();//new showTextureMapping();//new showProjectionPath(); //new showTensor();// new IS_Match_Compare();
	glutInit(&argc, argv);

	//glutInitContextFlags(GLUT_CORE_PROFILE);
	//glutInitContextVersion(4, 5);
	//cout << ">设置版本3.2" << endl;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(SCR_WIDTH, SCR_HEIGHT);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Tutorial 01");
	glutKeyboardFunc(keyCallback);
	
	GLenum res = glewInit();
	if (res != GLEW_OK)
	{
		printf("Error:'%s'\n", glewGetErrorString(res));
		return 1;
	}


	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	clock_t start_t, end_t;
	start_t = clock();
	currentScript->onInit();
	cout << "-------------------------------------------------------" << endl;
	end_t = clock();
	double time = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
	cout << "腳本初始化總計耗時:" << time << endl;
	//註冊鼠標funciton
	glutMouseFunc(&mouseCallBack);
	glutMotionFunc(&mouseDragCallBack);
	InitializeGlutCallbacks();
	//IMGUI_CHECKVERSION();
	//ImGui::CreateContext();
	//ImGuiIO& io = ImGui::GetIO(); (void)io;
	//ImGui::StyleColorsDark();
	//ImGui_ImplGLUT_Init();
	//ImGui_ImplGLUT_InstallFuncs();
	//ImGui_ImplOpenGL2_Init();
	glutMainLoop();
	//ImGui_ImplOpenGL2_Shutdown();
	//ImGui_ImplGLUT_Shutdown();
	//ImGui::DestroyContext();
	system("pause");
	return 0;
}

