# This file is meant to be used from within C++ and is thus written
# to be as easy-to-use as possible from the clunky C++ Python API.
# As much as possible is hard-coded.

import torch
import math
import numpy as np
import gymnasium as gym
from pathlib import Path

from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env, DummyVecEnv
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold

class SageRunner:
	def __init__(self):
		self.useIntAction = True
		self.model = []
		self.env = None
	
	# We want a simple 1-environment agent:
	def getEnv(self, renderMode=None):
		return DummyVecEnv([lambda: gym.make("CartPole-v1", render_mode=renderMode)])

	# Training is independent and won't be called from C++:
	def train(self, saveLocation):
		vec_env = self.getEnv()
		eval_env = gym.make("CartPole-v1")
	
		callback_on_best = StopTrainingOnRewardThreshold(reward_threshold=500, verbose=1)
		eval_callback = EvalCallback(eval_env, callback_on_new_best=callback_on_best, verbose=1)	

		model = PPO("MlpPolicy", vec_env, verbose=1)

		model.learn(total_timesteps=100000, callback=eval_callback)
		model.save(saveLocation)

	def loadMultiplePPOSages(self, directory : str):
		self.model = []
		searchDir = Path(directory)

		# Assume every zip file in the folder is a PPO model:
		for x in searchDir.iterdir():
			if str(x).endswith(".zip"):
				self.model.append(PPO.load(str(x)[:-4])) # Cut off the ".zip"

	def load(self):
		if Path("sageCommittee").is_dir():
			self.loadMultiplePPOSages("sageCommittee")
		else:
			self.model = [PPO.load('ppo_sage')]
	
		self.env = self.getEnv(renderMode='rgb_array')

	def render(self):		
		self.env.render()

	def getModelAction(self, obs):
		# Get the first action:
		action, _ = self.model[0].predict(obs)

		# Add in the rest (possibly no more):
		for i in range(1, len(self.model), 1):
			tmpAction, _ = self.model[i].predict(obs)
			action[0] += tmpAction[0]

		# Get the average as a vote:
		action[0] = int(round(action[0] / len(self.model), 0))
		return action

	def reset(self):
		# Gather observations and actions:
		obs = self.env.reset()
		actions = self.getModelAction(obs)
	
		# We only want the first observation and first action:
		actions = self.getList(actions[0])
		obs = self.getList(obs[0])
		
		return obs, actions

	# We often need to make sure our return value is a list. Do that here:
	def getList(self, inputVal):
		if isinstance(inputVal, np.ndarray):
			return list(inputVal)
		else:
			return list([inputVal])

	# Reminder, as a sage it wants to return an observation and what the sage-model
	# would do under those circumstances.
	def step(self, actualAction):
		if self.useIntAction:
			actualAction = int(actualAction)
		
		inAction = (actualAction, )
		obs, rewards, dones, info = self.env.step(inAction)
		actions = self.getModelAction(obs)
	
		obs = self.getList(obs[0])
		actions = self.getList(actions[0])
		rewards = rewards[0]
		dones = dones[0]
				
		# Vec-env returns need to be converted to what the singular
		# agent needs:
		return (obs, float(rewards), bool(dones), actions)

	def test(self):
		vec_env = self.getEnv(renderMode='human')
		self.load()

		obs = vec_env.reset()
		maxSteps = 10001
		# With a simple environment like this, what we want is in subscript
		# 0 of each of action, rewards, dones, and info.
		for i in range(1, maxSteps, 1):
			action = self.getModelAction(obs)
			obs, rewards, dones, info = vec_env.step(action)
			vec_env.render()
			if dones[0]:
				print(f"Step {i}: Done. Resetting the environment.")

def getInstance():
	return SageRunner()

if __name__ == "__main__":
	# runner = getInstance()
	
	# Train:
	# for i in range(7):
		# train(f"sageCommittee\\ppo_sage_{i}")
	
	# Test:
	# runner.test()
	
	print("In Python main.", flush=True)