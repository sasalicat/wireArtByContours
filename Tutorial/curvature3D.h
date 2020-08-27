#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
//用來計算曲率的類別,實際上沒有用到
class curvature3D
{
public:
	curvature3D();
	curvature3D(glm::vec3 p1,glm::vec3 p2,glm::vec3 p3);
	curvature3D(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3);
	~curvature3D();
	float getResult();
protected:
	glm::vec2 point2Ds[3];
};

