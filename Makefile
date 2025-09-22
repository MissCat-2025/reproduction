###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# RACCOON_DIR      - Root directory of the RACCOON project
#
###############################################################################

# 强制设置正确的路径，避免路径冲突
MOOSE_DIR := $(shell dirname `pwd`)/moose
RACCOON_DIR := $(shell dirname `pwd`)/raccoon

#
###############################################################################
# 注释掉MOOSE子模块检测，避免使用raccoon的moose子模块
# MOOSE_SUBMODULE    := $(CURDIR)/moose
# ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
#   MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
# else
#   MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
# endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).
# 排除特定的 MOOSE 源文件
ALL_MODULES                 := no

CHEMICAL_REACTIONS          := no
CONTACT                     := yes
ELECTROMAGNETICS            := no
EXTERNAL_PETSC_SOLVER       := no
FLUID_PROPERTIES            := no
FSI                         := no
FUNCTIONAL_EXPANSION_TOOLS  := no
GEOCHEMISTRY                := no
HEAT_TRANSFER               := yes
LEVEL_SET                   := no
MISC                        := no
NAVIER_STOKES               := no
OPTIMIZATION                := no
PERIDYNAMICS                := no
PHASE_FIELD                 := yes
POROUS_FLOW                 := no
RAY_TRACING                 := yes
REACTOR                     := no
RDG                         := no
RICHARDS                    := no
SOLID_MECHANICS             := yes
STOCHASTIC_TOOLS            := no
THERMAL_HYDRAULICS          := no
XFEM                        := no

include $(MOOSE_DIR)/modules/modules.mk

################################## RACCOON ####################################
# Include RACCOON as a dependency
RACCOON            := yes
RACCOON_LIB        := yes
# 添加RACCOON的所有include路径
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/materials/small_deformation_models
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/materials/hardening_models
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/interfaces
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/materials
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/kernels
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/base
ADDITIONAL_INCLUDES += -I$(RACCOON_DIR)/include/utils
ADDITIONAL_LIBS     += -L$(RACCOON_DIR)/lib -lraccoon-$(METHOD)


###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := reproduction
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
